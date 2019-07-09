#!/bin/bash
#SBATCH --job-name=extractReads
#SBATCH --output=sbatch.extractReads
#
#SBATCH --ntasks=40
#SBATCH --time=920:00
#SBATCH --mem=480000


if [ ! $# == 5 ]; then
  echo "Usage: $0 [Nanopore Reads BAM] [Chromium Reads] [GTF File] [Output Dir] [Merge factor]"
  exit
fi

mkdir $4

prefix=$4
prefix+="/reads"

date
executables/extractReads -b $2 -g $3 -o $prefix -t 40 -s $5 -su _chrom -e 3000 > extract.log
date
executables/extractReads -b $1 -g $3 -o $prefix -t 40 -s $5 -su _nano -ov >> extract.log
date
