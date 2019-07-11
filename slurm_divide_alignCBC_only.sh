#!/bin/bash
#SBATCH --job-name=a_CBC2
#SBATCH --output=sbatch.align2
#
#SBATCH --ntasks=40
#SBATCH --time=1440:00
#SBATCH --mem=100000

module load bedtools
module load samtools
module load flexbar

if [ ! $# == 5 ]; then
  echo "Usage: $0 [Path + Prefix of Nano and Chrom reads] [Output Prefix] [Start][End (number of Files)] [CellBarcodes]"
  exit
fi

leftT=$2
rightT=$2
P1_log=$2
leftT+="leftTail.log"
rightT+="rightTail.log"
P1_log+="P1.log"
# empty merge log file
echo "" > $leftT
echo "" > $rightT
#mkdir tmpfiles
tdir=$(mktemp -d tmp.XXXXXXXXX)
prefix=$1
st=$3
end=$4

declare -A barc
while read p; do
  barc["$p"]="1"
done <$5

for (( i=$st; i <= $end; i++)); do
  namechrom="${prefix}${i}_chrom.bam"
  namenano="${prefix}${i}_nano.bam"
  if [ ! -f $namechrom ] || [ ! -f $namenano ]; then
    continue
  fi
  echo "$namenano"
  nsize=$(samtools view $namenano | wc -l)
  csize=$(samtools view $namechrom | wc -l)
  echo "Number of Nanopore reads: "
  echo "$nsize"
  echo "Chromium reads: "
  echo "$csize"
  if [ $nsize == 0 ] || [ $csize == 0 ]; then
    echo "No reads nanopore or chromium reads found mapping selected exon."
    echo "Therefore skipping alignment."
    continue
  fi
  CBC_P=$namechrom
  CBC_P=${CBC_P%.*}
  CBC_="${CBC_P}_CBC.fasta"
  CBC_NAMES="${CBC_P}_CBC_names.fasta"
  CBC_TMP="${CBC_P}_CBC.tmp"
  CBC_TMP_CELL="${CBC_P}_CBC_CELL.tmp"


  NO_CDNA=$namenano
  NO_CDNA=${NO_CDNA%.*}
  NO_CDNA+="_no_cDNA.fasta"


  LEFT="${tdir}/leftTail.tmp"
  RIGHT="${tdir}/rightTail.tmp"
 
  # extract CBC from illumina reads
  samtools view $namechrom | perl -e 'while(<>){next unless /CB\:Z\:([ACGT]+)/; $CB=$1; /UB\:Z\:([ACGT]+)/; $UMI=$1; @tmp=split(/\t+/,$_); print ">",$tmp[0],"\n"; print $CB."\n";}' > $CBC_
  

  pp=$(wc -l "$CBC_" | cut -d " " -f 1)
  echo "CBC_ size: $pp"
  echo "Collect correct CBC"
  # 
  grep -v '>' $CBC_ |  awk '!a[$0]++'  > $CBC_TMP
  # | cut -b1-16|
  while read p; do
    if [[ "${barc["${p:0:16}"]}" == "1" ]];then
      echo "$p" >> $CBC_TMP_CELL
    fi
  done < $CBC_TMP

  date
  celll=$(wc -l "$CBC_TMP_CELL" | cut -d " " -f 1)
  echo "Selected barcodes from verified Cellbarcodes: "
  echo "$celll"

  #copy CBC into the name. Append unique number so flexbar has unique names for CBC combination
  perl -ne '$i++; $_ =~ s/[\r\n]+$//; print ">" . substr($_, 0, 16) . "|${i}\n" .  "$_\n"' $CBC_TMP_CELL > $CBC_NAMES
  #still need to remove empty lines due to very short sequences
  sed -i '/^$/d' $CBC_NAMES

  rm $CBC_TMP_CELL
  
  # remove cDna from nanoporereads using the alignment file
  executables/removecDNA -b $namenano -o $NO_CDNA

  # cut off P1 primer in both strand directions
  flexbar -r $NO_CDNA -t "${tdir}/P1" -as "CTACACGACGCTCTTCCGATCT" -n 40 -g -ae ANY -ac -ao 20 -l ALL -ag -1 -ai -2 -at 0.15

  cat ${tdir}/P1.log >> $P1_log 

  #split reads according to their direction
  executables/splitReads -f "${tdir}/P1.fasta" -o "${tdir}/P1_split" -l 25

  echo "$namechrom"
  date
  sizel=$(wc -l "${tdir}/P1_split_left_tail_trimmed.fasta" | cut -d " " -f 1)
  echo "P1_left $sizel"
  if [ $sizel -ne 0 ]; then
      executables/flexbar -r "${tdir}/P1_split_left_tail_trimmed.fasta" -t $LEFT -a $CBC_NAMES -n 40 -N 20 -g -ae LTAIL -ao 10 -l ALL -eve -ag -1 -ai -1 -at 0.24
  fi
  date
  sizer=$(wc -l ${tdir}/P1_split_right_tail_trimmed.fasta | cut -d " " -f 1)
  echo "P1right $sizer"
  if [ $sizer -ne 0 ]; then
      executables/flexbar -r "${tdir}/P1_split_right_tail_trimmed.fasta" -t $RIGHT -a $CBC_NAMES -n 40 -N 20 -g -ae RTAIL -ac -ao 10 -l ALL -eve -ag -1 -ai -1 -at 0.24
  fi
  date
#  echo "Copied Logs"
  cat "${LEFT}.log" >> $leftT
  cat "${RIGHT}.log" >> $rightT
  rm $CBC_
  rm $CBC_NAMES
  rm $CBC_TMP
  rm $NO_CDNA
done
rm $tdir -r
echo "Finished"
