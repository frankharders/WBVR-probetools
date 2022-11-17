#!/bin/bash

## design probe sequences 120bp long as default value.
## package is ProbeTools Stats v0.1.9 https://github.com/KevinKuchinski/ProbeTools
## several programs are install in the conda env
## https://stackabuse.com/how-to-parse-command-line-arguments-in-bash/


## date: 20221006-11:00
## add blat analysis to it
## add blat2gff.pl to it
## add rename the probe header from the probelist to a prefix with probetools conditions and numbers.
## In this way we exclude double probenames if probe output files are concatenated
## create probe hit table for metrics


##  activate the environment for this downstream analysis
eval "$(conda shell.bash hook)";
conda activate 02-viralprobes;


WORKDIR="$PWD";

mkdir -p "$PWD"/LOGS;

LOG="$PWD"/LOGS;

date

while getopts e:b::i::c:: flag
do
    case "${flag}" in
        e) extension=${OPTARG};;
        b) batch=${OPTARG};;
        i) identity=${OPTARG};;
		c) probecov=${OPTARG};;
		a) program=${OPTARG};;
    esac
done


[ -z "$extension" ] && echo "use -e for extension";
[ -z "$batch" ] && echo "use -b for batch size (default value: 100 probes per batch)";
[ -z "$identity" ] && echo "use -i for identity (default value: 95%)";
[ -z "$probecov" ] && echo "use -c for probe coverage";


find *.$extension > fasta.lst;




count0=1;
countS=$(cat fasta.lst | wc -l);

nodes=90;
lowcov=0;

for identity in 90 92 95 98 100;do

batch=50;
probecov=10;


while [ $count0 -le $countS ];do





	FILEin=$(cat fasta.lst | awk 'NR=='$count0 );

short=$(basename $FILEin .$extension); 

echo $short;
echo $FILEin;
OUTdir=$PWD/$short'-prbgroup'$batch'-cov'$probecov'-id'$identity;

rm -rf $OUTdir;
mkdir -p $OUTdir;

echo $OUTdir;

LOG1="$LOG"/"$short"."$batch"."$probecov"."$identity".makeprobes.log;
LOG2="$LOG"/"$short"."$batch"."$probecov"."$identity".getlowcov.log;
LOG3="$LOG"/"$short"."$batch"."$probecov"."$identity".stats.log;
LOG4="$LOG"/"$short"."$batch"."$probecov"."$identity".rename.log;

probetools makeprobes -t "$FILEin" -b "$batch" -o "$OUTdir"/"$short" -i "$identity" -T "$nodes" -D "$lowcov" -l 110 -c 100 -L 10 -d 5 > "$LOG1" 2>&1;

LOWCOVin="$OUTdir"/"$short"_capture.pt;

probetools getlowcov -i "$LOWCOVin" -o "$OUTdir"/"$short" > "$LOG2" 2>&1;


probetools stats -i "$LOWCOVin" -o "$OUTdir"/"$short" > "$LOG3" 2>&1;

RENAMEin="$OUTdir"/"$short"_probes.fa;
RENAMEout="$OUTdir"/"$short"_probes.renamed.fa;

rename.sh in="$RENAMEin" out="$RENAMEout" prefix="$short" ow > "$LOG4" 2>&1;

BLATdb="$FILEin";
BLATout="$OUTdir"/"$short".blat.out.psl;
GFFout="$OUTdir"/"$short".blat.out.gff;


PROBElst="$OUTdir"/"$short".probeName.lst;



		blat  "$BLATdb" "$RENAMEout" -t=dna -q=dna -noTrimA -out=psl -tileSize=11 -stepSize=1 -oneOff=2 -minIdentity=95 "$BLATout";

perl blat2gff.pl < "$BLATout" > "$GFFout";


cat "$RENAMEout" | grep '^>' | cut -f2 -d'>' > "$PROBElst";


PROBEhit="$OUTdir"/"$short".prb.hit.table.tab;

while read prb;do

tel=$(cat $BLATout | grep -c "$prb");
echo -e "$prb\t$tel" > "$PROBEhit";

done < "$PROBElst"


count0=$((count0+1));

done

done 

exit 1



### ProbeTools v0.1.9 https://github.com/KevinKuchinski/ProbeTools
### 
### Available modules:
### makeprobes - probe panel design using a general purpose incremental strategy
### clusterkmers - enumerate and cluster kmers from target sequences
### capture - assess probe panel coverage of target sequences
### getlowcov - extract low coverage sequences from target space
### stats - calculate probe coverage and depth stats
### merge - merge two sets of capture results


###   MODULE "makeprobes"
### Usage: probetools makeprobes -t <target seqs> -b <batch size> -o <path to output directory>/<output name> [<optional args>]
### Required arguments:
###  -t : path to target sequences in FASTA file
###  -b : number of probes in each batch (min=1)
###  -o : path to output directory and name to append to output files
### Optional arguments:
###  -m : max number of probes to add to panel (default=MAX, min=1)
###  -c : target for 10th percentile of probe coverage (default=90, min=1, max=100)
###  -k : length of probes to generate (default=120, min=32)
###  -s : number of bases separating each kmer (default=1, min=1)
###  -d : number of degenerate bases to permit in probes (default=0, min=0)
###  -i : nucleotide sequence identity (%) threshold used for kmer clustering and probe-target alignments (default=90, min=50, min=100)
###  -l : minimum length for probe-target alignments (default=60, min=1)
###  -D : minimum probe depth threshold used to define low coverage sub-sequences (default=0, min=0)
###  -L : minimum number of consecutive bases below probe depth threshold to define a low coverage sub-sequence (default=40, min=1)
###  -T : number of threads used by VSEARCH and BLASTn for clustering kmers and aligning probes to targets (default=MAX for VSEARCH, default=1 for BLASTn, min=1)

###   MODULE "getlowcov"
### Usage: probetools getlowcov -i <input file> -o <path to output directory>/<output name> [<optional args>]
### Required arguments:
###  -i : path to capture results in PT file
###  -o : path to output directory and name to append to output files
### Optional arguments:
###  -k : minimum sub-sequence length extracted, should be same as kmer length used for making probes (default=120, min=32)
###  -D : minimum probe depth threshold used to define low coverage sub-sequences (default=0, min=0)
###  -L : minimum number of consecutive bases below probe depth threshold to define a low coverage sub-sequence (default=40, min=1)


###   MODULE "stats"
### Usage: probetools stats -i <input file> -o <path to output directory>/<output name>
### Required arguments:
###  -i : path to capture results in PT file
###  -o : path to output directory and name to append to output files


###   BLAT analysis
#blat - Standalone BLAT v. 35 fast sequence search command line tool
#usage:
#   blat database query [-ooc=11.ooc] output.psl
#where:
#   database and query are each either a .fa , .nib or .2bit file,
#   or a list these files one file name per line.
#   -ooc=11.ooc tells the program to load over-occurring 11-mers from
#               and external file.  This will increase the speed
#               by a factor of 40 in many cases, but is not required
#   output.psl is where to put the output.
#   Subranges of nib and .2bit files may specified using the syntax:
#      /path/file.nib:seqid:start-end
#   or
#      /path/file.2bit:seqid:start-end
#   or
#      /path/file.nib:start-end
#   With the second form, a sequence id of file:start-end will be used.
#options:
#   -t=type     Database type.  Type is one of:
#                 dna - DNA sequence
#                 prot - protein sequence
#                 dnax - DNA sequence translated in six frames to protein
#               The default is dna
#   -q=type     Query type.  Type is one of:
#                 dna - DNA sequence
#                 rna - RNA sequence
#                 prot - protein sequence
#                 dnax - DNA sequence translated in six frames to protein
#                 rnax - DNA sequence translated in three frames to protein
#               The default is dna
#   -prot       Synonymous with -t=prot -q=prot
#   -ooc=N.ooc  Use overused tile file N.ooc.  N should correspond to
#               the tileSize
#   -tileSize=N sets the size of match that triggers an alignment.
#               Usually between 8 and 12
#               Default is 11 for DNA and 5 for protein.
#   -stepSize=N spacing between tiles. Default is tileSize.
#   -oneOff=N   If set to 1 this allows one mismatch in tile and still
#               triggers an alignments.  Default is 0.
#   -minMatch=N sets the number of tile matches.  Usually set from 2 to 4
#               Default is 2 for nucleotide, 1 for protein.
#   -minScore=N sets minimum score.  This is the matches minus the
#               mismatches minus some sort of gap penalty.  Default is 30
#   -minIdentity=N Sets minimum sequence identity (in percent).  Default is
#               90 for nucleotide searches, 25 for protein or translated
#               protein searches.
#   -maxGap=N   sets the size of maximum gap between tiles in a clump.  Usually
#               set from 0 to 3.  Default is 2. Only relevent for minMatch > 1.
#   -noHead     suppress .psl header (so it's just a tab-separated file)
#   -makeOoc=N.ooc Make overused tile file. Target needs to be complete genome.
#   -repMatch=N sets the number of repetitions of a tile allowed before
#               it is marked as overused.  Typically this is 256 for tileSize
#               12, 1024 for tile size 11, 4096 for tile size 10.
#               Default is 1024.  Typically only comes into play with makeOoc.
#               Also affected by stepSize. When stepSize is halved repMatch is
#               doubled to compensate.
#   -mask=type  Mask out repeats.  Alignments won't be started in masked region
#               but may extend through it in nucleotide searches.  Masked areas
#               are ignored entirely in protein or translated searches. Types are
#                 lower - mask out lower cased sequence
#                 upper - mask out upper cased sequence
#                 out   - mask according to database.out RepeatMasker .out file
#                 file.out - mask database according to RepeatMasker file.out
#   -qMask=type Mask out repeats in query sequence.  Similar to -mask above but
#               for query rather than target sequence.
#   -repeats=type Type is same as mask types above.  Repeat bases will not be
#               masked in any way, but matches in repeat areas will be reported
#               separately from matches in other areas in the psl output.
#   -minRepDivergence=NN - minimum percent divergence of repeats to allow
#               them to be unmasked.  Default is 15.  Only relevant for
#               masking using RepeatMasker .out files.
#   -dots=N     Output dot every N sequences to show program's progress
#   -trimT      Trim leading poly-T
#   -noTrimA    Don't trim trailing poly-A
#   -trimHardA  Remove poly-A tail from qSize as well as alignments in
#               psl output
#   -fastMap    Run for fast DNA/DNA remapping - not allowing introns,
#               requiring high %ID. Query sizes must not exceed 5000.
#   -out=type   Controls output file format.  Type is one of:
#                   psl - Default.  Tab separated format, no sequence
#                   pslx - Tab separated format with sequence
#                   axt - blastz-associated axt format
#                   maf - multiz-associated maf format
#                   sim4 - similar to sim4 format
#                   wublast - similar to wublast format
#                   blast - similar to NCBI blast format
#                   blast8- NCBI blast tabular format
#                   blast9 - NCBI blast tabular format with comments
#   -fine       For high quality mRNAs look harder for small initial and
#               terminal exons.  Not recommended for ESTs
#   -maxIntron=N  Sets maximum intron size. Default is 750000
#   -extendThroughN - Allows extension of alignment through large blocks of N's




