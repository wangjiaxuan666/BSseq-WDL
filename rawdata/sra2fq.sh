#!/bin/sh

func() {
    echo "Usage:"
    echo "[-i]: the input file .sra"
    echo "[-h]: The help document"
    echo "[-o]: the output file .fastq"
    exit -1
}

while getopts "i:o:h" opt; do
    case $opt in
      i) input="$OPTARG";;
      o) output="$OPTARG";;
      h) func;;
      ?) func;;
    esac
done

fasterq-dump -e 2 --split-files -O ./ --outfile $output $input
