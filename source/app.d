import bio.bam.pileup;
import bio.bam.reader;
import std.algorithm: filter, map, sort, uniq;
import std.array: split;
import std.conv: to;
import std.file: exists;
import std.format;
import std.getopt;
import std.range;
import std.stdio;

static string usage =
"ac - alleleCounter clone\n"~
"------------------------\n\n"~
"Scans a BAM alignment for a number of genomic coordinates and reports the number of As, Cs, Gs and Ts at each.\n\n"~
"Basic Usage:\n"~
"    ac -b|--bamfile <bamfile> -l|--locifile <locifile>\n"~
"    ac -h|--help\n\n"~
"Available options:\n"~
"    -b|--bamfile       - The BAM alignment file to search\n"~
"    -l|--locifile      - A tab-separated list of Chromosome-Position pairs\n"~
"    -m|--minmapqual    - Reads with mapping quality below this value are skipped (default: 35)\n"~
"    -q|--minbasequal   - Bases below this quality value are ignored (default: 20)\n"~
"    -f|--required-flag - Reads must match all flags contained in this integer (default: 3)\n"~
"    -F|--filtered-flag - Reads must not match any flags contained in this integer (default: 3852)\n";

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
    int fflags = 3852;
    int rflags = 3;

    try {
        auto args = getopt(
                argv,
                std.getopt.config.required, "bamfile|b", "Path to sample BAM file.", &bamfile,
                std.getopt.config.required, "locifile|l", "Path to loci file.", &locifile,
                "minbasequal|m", "Minimum base quality [Default: 20].", &minbasequal,
                "minmapqual|q", "Minimum mapping quality [Default: 35].", &minmapqual,
                config.caseSensitive, "filtered-flag|F", "Ignore reads matching flags [Default: 3852].", &fflags,
                config.caseSensitive, "required-flag|f", "Only count reads matching flags [Default: 3].", &rflags);

        if (args.helpWanted) {
            defaultGetoptPrinter(usage, args.options);
            return;
        }
    }
    catch (GetOptException) {
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

    writefln("#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth");

    foreach (linenum, line; loci.byLineCopy.enumerate(1)) {
        auto spl = split(line, '\t');
        string refname = to!string(spl[0]);
        ulong pos_1based;

        uint ref_id;
        try {
            ref_id = bam.ref_id(refname);
        }
        catch (Throwable) {
            stderr.writefln("Skipping locus on line %d with unknown chromosome: %s", linenum, line);
            continue;
        }

        try {
            pos_1based = to!ulong(spl[1]);
        }
        catch (std.conv.ConvException) {
            stderr.writefln("Skipping locus on line %d with undecipherable position: %s", linenum, line);
            continue;
        }

        auto pos_0based = pos_1based - 1;

        if (ref_id != curr_ref) {
            curr_ref = ref_id;
            pileup = makePileup(bam.reference(curr_ref)[1 .. uint.max]);
            column = pileup.front;
        }

        if (pileup.empty) {
            writefln("%s\t%d\t0\t0\t0\t0\t0", refname, pos_1based);
            continue;
        }

        assert(column.ref_id == ref_id);
        while(column.position < pos_0based && column.ref_id == ref_id) {
            if (pileup.empty) {
                writefln("%s\t%d\t0\t0\t0\t0\t0", refname, pos_1based);
                break;
            }
            pileup.popFront();
            column = pileup.front;
        }

        if (column.position == pos_0based) {
            auto bases = column.reads
                .filter!(read => (read.current_base_quality >= minbasequal) && (read.mapping_quality >= minmapqual) && (read.flag() & fflags) == 0 && (read.flag() & rflags) == rflags)
                .array
                .sort!"a.name < b.name"
                .uniq!"a.name == b.name"
                .map!(read => read.current_base)
                .to!string;
            writefln("%s\t%d\t%s", refname, pos_1based, count_bases(bases));
        }

        if (column.position > pos_0based) {
            writefln("%s\t%d\t0\t0\t0\t0\t0", refname, pos_1based);
            continue;
        }
    }
}

