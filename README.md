# ac
Basic alleleCounter implementation

Build:

    dub build --compiler=ldc2 --build=release
    
Run:

    ./ac -b data.bam -l loci.txt  --minbasequal 20 --minmapqual 35
    
