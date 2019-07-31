#!/bin/bash
#SBATCH --job-name=a_CBC
#SBATCH --output=sbatch.align
#
#SBATCH --ntasks=40
#SBATCH --time=1440:00
#SBATCH --mem=100000

module load samtools
module load flexbar/3.5.0

if [ ! $# == 3 ]; then
  echo "Usage: $0 [NanoporeReads bam][CellBarcodes Whitelist] [Output Prefix]"
  exit
fi

leftT=$3
rightT=$3
P1_log=$3
leftT+="leftTail.log"
rightT+="rightTail.log"
P1_log+="P1.log"
# empty merge log file
echo "" > $leftT
echo "" > $rightT
#mkdir tmpfiles
tdir=$(mktemp -d tmp.XXXXXXXXX)
namenano=$1
whitelist=$2

CBC="${tdir}/whitelist.fasta"

NO_CDNA="${tdir}/nano_no_cDNA.fasta"

LEFT="${tdir}/leftTail.tmp"
RIGHT="${tdir}/rightTail.tmp"
 
#copy CBC into the name. Append unique number so flexbar has unique names for CBC combination
perl -ne '$i++; $_ =~ s/[\r\n]+$//; print ">" . substr($_, 0, 16) . "\n" .  "$_\n"' $whitelist > $CBC
#still need to remove empty lines due to very short sequences
#sed -i '/^$/d' $CBC_NAMES

  
# remove cDna from nanoporereads using the alignment file
/scratch/sboenigk/executables/single_cell_scripts_exe/removecDNA -b $namenano -o $NO_CDNA

# cut off P1 primer in both strand directions
/scratch/sboenigk/executables/lastest_flexbar/flexbar -r $NO_CDNA -t "${tdir}/P1" -as "CTACACGACGCTCTTCCGATCT" -n 40 -g -at ANY -ac ON -ao 20 -l ALL -ag -1 -ai -2 -ae 0.15

cat ${tdir}/P1.log >> $P1_log 


#split reads according to their direction
/scratch/sboenigk/executables/single_cell_scripts_exe/splitReads -f "${tdir}/P1.fasta" -o "${tdir}/P1_split" -c -l 25
date
sizel=$(wc -l "${tdir}/P1_split_left_tail_trimmed.fasta" | cut -d " " -f 1)
echo "P1_left $sizel"
if [ $sizel -ne 0 ]; then
  /scratch/sboenigk/executables/lastest_flexbar/flexbar -r "${tdir}/P1_split_left_tail_trimmed.fasta" -t $LEFT -a $CBC -n 40 -N 20 -g -at LTAIL -ao 10 -l ALL -alt -ag -1 -ai -1 -ae 0.24
fi
date
sizer=$(wc -l ${tdir}/P1_split_right_tail_trimmed.fasta | cut -d " " -f 1)
echo "P1right $sizer"
if [ $sizer -ne 0 ]; then
  /scratch/sboenigk/executables/lastest_flexbar/flexbar -r "${tdir}/P1_split_right_tail_trimmed.fasta" -t $RIGHT -a $CBC -n 40 -N 20 -g -at RTAIL -ac ONLY -ao 10 -l ALL -alt -ag -1 -ai -1 -ae 0.24
fi
date
#  echo "Copied Logs"
cat "${LEFT}.log" >> $leftT
cat "${RIGHT}.log" >> $rightT
#rm $CBC
#rm $NO_CDNA
#rm $tdir -r
echo "Finished"
