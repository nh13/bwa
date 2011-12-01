FQ="ce_small.fq.gz"
FA="Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.gz"
BWA="../bwa"

usage() {
	echo -e "\n\nWarning:$1"
	echo -e "\nUSAGE:\n$0 [-a path/to/alernative/bwa] [seq.fq = $FQ] \\\\"
	echo -e "\t[ref.fa = $FA]\n\t(without -a $BWA is used)\n"
	[ -f ce_med.fq.gz ] echo -e "use ce_med.fq.gz as fastq to test more reads..."
	exit;
}

[ $# -gt 0 -a "$1" = "-a" ] && {
	TAG="_alt";
	shift;
	[ -x "$1" ] || { usage "$1 is not executable"; exit; }
	BWA="$1"
	shift;
}
[ $# -gt 0 ] && FQ="$2"
[ $# -gt 1 ] && FA="$3";

dest="./out"
src="."
echo using $BWA

BN="`basename $FQ ".fq.gz"`"

CWD="`pwd`"
cd ".."
[ -f Makefile -a -f bwa -a `date --reference=Makefile +%s%N` -lt \
  `date --reference=bwa +%s%N` ] || make clean
make
[ $? -ne 0 ] && { usage "make failed:$!"; }
cd "$CWD"
echo "----------aln"
$BWA aln -t 4 $src/$FA ${FQ} > $dest/${BN}${TAG}.sai

echo "---------- sampe"
$BWA samse $src/$FA $dest/${BN}${TAG}.sai ${FQ} > $dest/${BN}${TAG}.sam

echo "testing output ..."
#ignore first characters, contains bwa version
[ -f $dest/${BN}_alt.sam -a -f $dest/${BN}.sam ] || { usage "$dest/${BN}(_alt).sam does not yet exist"; }
cmp -i2850 $dest/${BN}.sam $dest/${BN}_alt.sam
echo "(silence is good in this case)"
