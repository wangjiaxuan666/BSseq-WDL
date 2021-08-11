#!/bin/sh
func() {
    echo "Usage:"
    echo "[-m]: the manifest which have the sampleid,read1,read2,the file is \t sep table"
    echo "[-h]: The help document"
    exit -1
}

manifest="./manifest"

while getopts "m:h" opt; do
    case $opt in
      m) manifest="$OPTARG";;
      h) func;;
      ?) func;;
    esac
done

cat $manifest | while read col1 col2 col3
do
  sampleid=$col1
  read1=$col2
  read2=$col3
  mkdir $sampleid
  cd $sampleid
  fastp -i $read1 \
    -I $read2 \
    -o ./$sampleid.filter.R1.fq.gz \
    -O ./$sampleid.filter.R2.fq.gz \
    -R './$sampleid_fastp_QC_report'
  read1path=$(ls *.filter.R1.fq.gz| sed "s:^:`pwd`/:")
  read2path=$(ls *.filter.R2.fq.gz| sed "s:^:`pwd`/:")
  rm ../afterqc.manifest
  echo -e "$sampleid\t$read1path\t$read2path" >> ../afterqc.manifest
  cd ../
  echo "=================="
  echo "$sampleid have fastp QC filtered"
  echo "=================="
done

