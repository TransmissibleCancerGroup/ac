import bio.bam.pileup;
import bio.bam.reader;
import std.algorithm: filter, map, sort, chunkBy;
import std.array: split;
import std.conv: to, ConvException;
import std.file: exists;
import std.format;
import std.getopt;
import std.range;
import std.stdio;

static string usage = "ac - alleleCounter clone\nUsage: ac -b|--bamfile <bamfile> -l|--locifile <locifile> -u|--umi\n       ac -h|--help";

int ref_id(R)(R reader, string ref_name) {
    return reader[ref_name].id;
}

auto count_bases(string bases) {
    int nA, nC, nG, nT;
    foreach(base; bases) {
        switch(base) {
            case 'A': nA++; break;
            case 'C': nC++; break;
            case 'G': nG++; break;
            case 'T': nT++; break;
            default: break;
        }
    }
    auto s = format!"%d\t%d\t%d\t%d\t%d"(nA, nC, nG, nT, nA+nC+nG+nT);
    return s;
}

auto qual_to_ascii(Array)(Array qualities) {
    map!(v => v+33).array();
}

void main(string[] argv)
{
    string bamfile;
    string locifile;
    int minmapqual = 35;
    int minbasequal = 20;
    bool umi = false;

    try {
        auto args = getopt(
                argv,
                std.getopt.config.required, "bamfile|b", "Path to sample BAM file.", &bamfile,
                std.getopt.config.required, "locifile|l", "Path to loci file.", &locifile,
                "minbasequal|m", "Minimum base quality [Default: 20].", &minbasequal,
                "minmapqual|q", "Minimum mapping quality [Default: 35].", &minmapqual,
                "umi|u", "Output counts on per-UMI basis [Default: false]", &umi);

        if (args.helpWanted) {
            defaultGetoptPrinter(usage, args.options);
            return;
        }
    }
    catch (GetOptException) {
        writeln(usage);
        return;
    }
    catch (ConvException e) {
        writefln("Error understanding command line arguments: \"%s\"", e.msg);
        writeln(usage);
        return;
    }

    if(!exists(bamfile)) {
        writefln("File %s does not exist: exiting.", bamfile);
        return;
    }

    if(!exists(locifile)) {
        writefln("File %s does not exist: exiting.", locifile);
        return;
    }

    auto bam = new BamReader(bamfile);
    auto loci = File(locifile);
    scope(exit) {
        loci.close();
    }

    if (!bam.has_index()) {
        bam.createIndex();
    }

    int curr_ref = 0;
    auto pileup = makePileup(bam.reference(curr_ref)[1 .. uint.max]);
    auto column = pileup.front;

    if (umi) {
        writefln("#CHR\tPOS\tBarcode\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth\tStrand");
    }
    else {
        writefln("#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth");
    }

    foreach (line; loci.byLineCopy) {
        auto spl = split(line, '\t');
        string refname = to!string(spl[0]);
        auto ref_id = bam.ref_id(refname);
        ulong pos_1based = to!ulong(spl[1]);
        auto pos_0based = pos_1based - 1;

        if (ref_id != curr_ref) {
            curr_ref = ref_id;
            pileup = makePileup(bam.reference(curr_ref)[1 .. uint.max]);
            column = pileup.front;
        }

        if (pileup.empty) {
            if (umi) {
                writefln("%s\t%d\t.\t0\t0\t0\t0\t0\t.", refname, pos_1based);
            }
            else {
                writefln("%s\t%d\t0\t0\t0\t0\t0", refname, pos_1based);
            }
            continue;
        }

        assert(column.ref_id == ref_id);
        while(column.position < pos_0based && column.ref_id == ref_id) {
            if (pileup.empty) {
                if (umi) {
                    writefln("%s\t%d\t.\t0\t0\t0\t0\t0\t.", refname, pos_1based);
                }
                else {
                    writefln("%s\t%d\t0\t0\t0\t0\t0", refname, pos_1based);
                }
                break;
            }
            pileup.popFront();
            column = pileup.front;
        }

        if (column.position == pos_0based) {
            if (umi) {
                    auto forward_chunks = column.reads
                    .sort!((a, b) => a.sequence.to!string[$-5 .. $] < b.sequence.to!string[$-5 .. $])
                    .filter!(read => !read.is_reverse_strand())
                    .chunkBy!((a, b) => a.sequence.to!string[$-5 .. $] == b.sequence.to!string[$-5 .. $]);

                auto reverse_chunks = column.reads
                    .sort!((a, b) => a.sequence.to!string[0 .. 5] < b.sequence.to!string[0 .. 5])
                    .filter!(read => read.is_reverse_strand())
                    .chunkBy!((a, b) => a.sequence.to!string[0 .. 5] == b.sequence.to!string[0 .. 5]);
                
                foreach (chunk; forward_chunks) {
                    auto arr = chunk.array;
                    auto bases = arr
                        .filter!(read => (read.current_base_quality >= minbasequal) && (read.mapping_quality >= minmapqual))
                        .map!(read => read.current_base)
                        .to!string;
                    writefln("%s\t%d\t%s\t%s\t+", refname, pos_1based, arr[0].sequence[$-5..$], count_bases(bases));
                }

                foreach (chunk; reverse_chunks) {
                    auto arr = chunk.array;
                    auto bases = arr
                        .filter!(read => (read.current_base_quality >= minbasequal) && (read.mapping_quality >= minmapqual))
                        .map!(read => read.current_base)
                        .to!string;
                    writefln("%s\t%d\t%s\t%s\t-", refname, pos_1based, arr[0].sequence[0..5], count_bases(bases));
                }
            // TODO: if reverse(UMI) of reverse strand == UMI of forward strand, collapse both outputs onto
            //       a single line by adding up the bases
            // TODO: tool should filter out UMI duplicates, so only report the first of each UMI set and skip the rest,
            //       or even collapse everything down so there is one line per position
            }
            else {
                auto bases = column.reads
                    .filter!(read => (read.current_base_quality >= minbasequal) && (read.mapping_quality >= minmapqual))
                    .map!(read => read.current_base)
                    .to!string;
                writefln("%s\t%d\t%s", refname, pos_1based, count_bases(bases));
            }
        }

        if (column.position > pos_0based) {
            if (umi) {
                writefln("%s\t%d\t.\t0\t0\t0\t0\t0\t.", refname, pos_1based);
            }
            else {
                writefln("%s\t%d\t0\t0\t0\t0\t0", refname, pos_1based);
            }
            continue;
        }
    }
}

