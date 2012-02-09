FQ="ce_small.fq.gz"
FA="Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.gz"
DIR=".."

usage() {
	echo -e "\n\nWarning:$1"
	echo -e "\nUSAGE:\n$0 [-a path/to/alernative/bwa] [seq.fq = $FQ] \\\\"
	echo -e "\t[ref.fa = $FA]\n\t(without -a $BWA is used)\n"
	[ -f ce_med.fq.gz ] && echo -e "use ce_med.fq.gz as fastq to test more reads..."
	exit;
}

[ $# -gt 0 -a "$1" = "-a" ] && {
	TAG="_alt";
	shift;
	DIR="$1"
	[ -d "$DIR" ] || DIR="${DIR/\/bwa\/bwa/\/bwa}"
	shift;
}
[ $# -gt 0 ] && FQ="$1"
[ $# -gt 1 ] && FA="$2";
[ -d $DIR ] || { usage "specify bwa directory with -a"; }
BWA="$DIR/bwa"

CWD="`pwd`"
dest="$CWD/out"
log="$CWD/log"
[ -z "`which gprof`" ] && { usage "gprof is not installed"; }
echo using $BWA

BN="`basename $FQ ".fq.gz"`"

cd "$DIR"
patch -p1 < "$CWD/gprof_Makefile.patch"

make clean; make
patch -R -p1 < "$CWD/gprof_Makefile.patch"
[ $? -ne 0 ] && { usage "make failed:$!"; }
cd "$CWD"

find_newf() {
	local i=0
        while [ -f "$1${i}$2" ]; do
		i=$((i+1))
	done
	echo ${i};
}

profile() {
	mv gmon.out "$DIR/"
	cd "$DIR"
	cp gmon.out $log/gmon/${1}_${BN}${TAG}_${2}_gmon.out
	gprof ./bwa > "$log/${1}_${BN}${TAG}_${2}_gprof.txt"
	for f in *.c *.h; do
		gcov -b $f
		[ -f "$f.gcov" ] && mv "$f.gcov" "$log/gcov/${f}_${1}_${BN}${TAG}_${2}.gcov"
	done
	cd $CWD
}

NDX=$(find_newf "$log/aln_${BN}${TAG}_" "_gprof.txt")

echo "----------aln"
$BWA aln -t 4 $CWD/$FA ${FQ} > $dest/${BN}${TAG}.sai
profile aln${TAG} $NDX

echo "---------- sampe"
$BWA samse $CWD/$FA $dest/${BN}${TAG}.sai ${FQ} > $dest/${BN}${TAG}.sam
profile samse${TAG} $NDX

echo "testing output ..."
#ignore first characters, contains bwa version
[ -f $dest/${BN}_alt.sam -a -f $dest/${BN}.sam ] || { usage "$dest/${BN}(_alt).sam does not yet exist"; }
cmp -i2850 $dest/${BN}.sam $dest/${BN}_alt.sam
echo "(silence is good in this case)"
#rm $dest/${BN}${TAG}{.sam,{1,2}.sai}
#cd "$DIR"
#rm *.gcda *.gcov *.gcno gmon.out
echo "output is in $log/ and $log/gcov/"

