#!/bin/bash
#SBATCH --job-name=a_CBC_UMI3
#SBATCH --output=sbatch.align3
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

#echo "Log: " > log.txt
# empty merge log file
leftT=$2
rightT=$2
P1_log=$2
leftT+="leftTail.log"
rightT+="rightTail.log"
P1_log+="P1.log"
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
  CBC_UMIP=$namechrom
  CBC_UMIP=${CBC_UMIP%.*}
  CBC_UMI="${CBC_UMIP}_CBC_UMI.fasta"
  CBC_UMI_NAMES="${CBC_UMIP}_CBC_UMI_names.fasta"
  CBC_UMI_TMP="${CBC_UMIP}_CBC_UMI.tmp"
  CBC_UMI_TMP_CELL="${CBC_UMIP}_CBC_UMI_CELL.tmp"


  NO_CDNA=$namenano
  NO_CDNA=${NO_CDNA%.*}
  NO_CDNA+="_no_cDNA.fasta"


  LEFT="${tdir}/leftTail.tmp"
  RIGHT="${tdir}/rightTail.tmp"

  CBC_only=false 
 
  # extract CBC + UMI from chromium reads
  samtools view $namechrom | perl -e 'while(<>){next unless /CB\:Z\:([ACGT]+)/; $CB=$1; /UB\:Z\:([ACGT]+)/; $UMI=$1; @tmp=split(/\t+/,$_); print ">",$tmp[0],"\n"; print $CB.$UMI."\n";}' > $CBC_UMI
  

  pp=$(wc -l "$CBC_UMI" | cut -d " " -f 1)
  echo "CBC_UMI size: $pp"
  if [ $pp -gt 30000 ]; then
    echo "Use CBC white list"
    cp $5 $CBC_UMI_TMP_CELL
    CBC_only=true
  else
    echo "Collect correct CBC"  
    #copy CBC into the name. Append unique number so flexbar has unique names for CBC + UMI combination
    grep -v '>' $CBC_UMI | sort -u > $CBC_UMI_TMP
    # | cut -b1-16|
    while read p; do
#    grep "^${p}" $CBC_UMI_TMP >> $CBC_UMI_TMP_CELL
      if [[ "${barc["${p:0:16}"]}" == "1" ]];then
        echo "$p" >> $CBC_UMI_TMP_CELL
      fi
    done < $CBC_UMI_TMP

    date
    celll=$(wc -l "$CBC_UMI_TMP_CELL" | cut -d " " -f 1)
    echo "Selected barcodes from verified Cellbarcodes: "
    echo "$celll"

    if [ $celll -gt 15000 ]; then
      CBC_only=true
    fi
  fi
  cp $CBC_UMI_TMP_CELL $CBC_UMI_TMP
  if [[ "$CBC_only" == "true" ]];then
      echo "CBC only"
     cat $CBC_UMI_TMP_CELL | cut -b1-16 | sort -u > $CBC_UMI_TMP 
  fi 

  perl -ne '$i++; $_ =~ s/[\r\n]+$//; print ">" . substr($_, 0, 16) . "|${i}\n" .  "$_\n"' $CBC_UMI_TMP > $CBC_UMI_NAMES
  #still need to remove empty lines due to very short sequences
  sed -i '/^$/d' $CBC_UMI_NAMES

  rm $CBC_UMI_TMP_CELL
  
#  echo "extracted CBC UMI"
  # remove cDna from nanoporereads using the alignment file
  executables/removecDNA -b $namenano -o $NO_CDNA

#  sizenc=$(wc -l "$NO_CDNA" | cut -d " " -f 1)
#  echo "nocDna $sizenc"
  
#  echo "nocDNA"

  # cut off P1 primer in both strand directions
  flexbar -r $NO_CDNA -t "${tdir}/P1" -as "CTACACGACGCTCTTCCGATCT" -g -ae ANY -ac -ao 20 -l ALL -ag -1 -ai -2 -at 0.2

  cat ${tdir}/P1.log >> $P1_log 

  #split reads according to their direction
  if [[ "$CBC_only" == "true" ]]; then
    executables/splitReads -f "${tdir}/P1.fasta" -o "${tdir}/P1_split" -l 25
  else
    echo "CBC + UMI extracted"
    executables/splitReads -f "${tdir}/P1.fasta" -o "${tdir}/P1_split" -l 35
  fi

  echo "$namechrom"
  date
  sizel=$(wc -l "${tdir}/P1_split_left_tail_trimmed.fasta" | cut -d " " -f 1)
  echo "P1_left $sizel"
  if [ $sizel -ne 0 ]; then
    if [[ "$CBC_only" == "true" ]]; then
      executables/lastest_flexbar/flexbar -r "${tdir}/P1_split_left_tail_trimmed.fasta" -t $LEFT -a $CBC_UMI_NAMES -n 40 -N 20 -g -ae LTAIL -ao 10 -l ALL -eve -ag -1 -ai -1 -at 0.24
    else
      executables/lastest_flexbar/flexbar -r "${tdir}/P1_split_left_tail_trimmed.fasta" -t $LEFT -a $CBC_UMI_NAMES -n 40 -N 2 -g -ae LTAIL -ao 20 -l ALL -eve -ag -1 -ai -1 -at 0.24
    fi
  fi
  date
  sizer=$(wc -l ${tdir}/P1_split_right_tail_trimmed.fasta | cut -d " " -f 1)
  echo "P1right $sizer"
  if [ $sizer -ne 0 ]; then
    if [[ "$CBC_only" == "true" ]]; then
      executables/lastest_flexbar/flexbar -r "${tdir}/P1_split_right_tail_trimmed.fasta" -t $RIGHT -a $CBC_UMI_NAMES -n 40 -N 20 -g -ae RTAIL -ac -ao 10 -l ALL -eve -ag -1 -ai -1 -at 0.24
    else
      executables/lastest_flexbar/flexbar -r "${tdir}/P1_split_right_tail_trimmed.fasta" -t $RIGHT -a $CBC_UMI_NAMES -n 40 -N 2 -g -ae RTAIL -ac -ao 20 -l ALL -eve -ag -1 -ai -1 -at 0.24
    fi
  fi
  date

#  echo "Copied Logs"
  cat "${LEFT}.log" >> $leftT
  cat "${RIGHT}.log" >> $rightT
  rm $CBC_UMI
  rm $CBC_UMI_NAMES
  rm $CBC_UMI_TMP
  rm $NO_CDNA
done
rm $tdir -r
echo "Finished"
