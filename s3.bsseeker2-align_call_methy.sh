#!/bin/sh

func() {
    echo "Usage:"
    echo "[-m]: the manifest which have the sampleid,read1,read2,the file is \t sep table"
    echo "[-h]: The help document"
    echo "[-g]: The split fq reads number the default is 10000000(10M)"
    echo "[-i]: The index dirpath which bs_seeker2-build teh reference genome index ,default is bowtie2"
    exit -1
}

manifest="./afterqc.manifest"
index_position="/hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/BSseeker/bs_utils/reference_genomes/genome.fa_bowtie2"
splitnum="20000000"

while getopts "i:m:g:h" opt; do
    case $opt in
      i) index_position="$OPTARG";;
      m) manifest="$OPTARG";;
      g) splitnum="$OPTARG";;
      h) func;;
      ?) func;;
    esac
done

#==============================================
# if first analysis need check the genome index
#==============================================

#/hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/miniconda3/envs/py27/bin/python \
#/hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/BSseeker/bs_seeker2-build.py \
#-f ./refer/genome.fa \
#--aligner=bowtie2

# conda activate env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate py27

cat $manifest | while read col1 col2 col3
do 
  sampleid=$col1
  read1=$col2
  read2=$col3
  if [ ! -d $sampleid ]; then mkdir $sampleid; fi
  cd $sampleid
  mkdir bam # save the bsseeker mapping bam files
  mkdir align_result # save the call mythy site infomation files
  mkdir log # save the log files

  zcat $read1 | split -l $splitnum - R1_
  zcat $read2 | split -l $splitnum - R2_
  #forloop2
  ls R1_[a-z][a-z]* | while read id
  do
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/miniconda3/envs/py27/bin/python \
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/BSseeker/bs_seeker2-align.py \
    -i $id \
    -g ../refer/genome.fa \
    --aligner=bowtie2 \
    --bt2-p 4 --bt2--end-to-end --bt2--very-sensitive --bt2--dovetail \
    --temp_dir=`pwd` \
    -m 0.1 \
    --XSteve \
    -o $id.bam
  done

  echo -e "\nthe $sampleid read1 mapping work have done!"

  ls R2_[a-z][a-z]* | while read id
  do
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/miniconda3/envs/py27/bin/python \
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/BSseeker/Antisense.py \
    -i $id \
    -o $id.anti.fq

    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/miniconda3/envs/py27/bin/python \
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/BSseeker/bs_seeker2-align.py \
    -i $id.anti.fq \
    -g ../refer/genome.fa \
    --aligner=bowtie2 \
    --bt2-p 8 --bt2--end-to-end --bt2--very-sensitive --bt2--dovetail \
    --temp_dir=`pwd` \
    -m 0.1 \
    --XSteve \
    -o $id.bam
  done
  
  echo "\nthe $sampleid read2 mapping work have done!"

  samtools merge merge.bam `ls R[1-2]_[a-z][a-z]*.bam`

  echo "\nthe $sampleid mapping bam files have merged!"

  if [ -a merge.bam ]
  then
      rm R[1-2]_[a-z][a-z]*.bam
      mv *.bam.bs_seeker2_log ./log
      rm R[1-2]_[a-z][a-z]*
  fi

  samtools sort merge.bam -o ./bam/sort.bam
  # 标签XS:i:1代表BS-Seeker2利用CH位点的甲基化信息判断出reads failed in bisulfite-conversion
  samtools view -h ./bam/sort.bam | grep -v 'XS:i:1' | samtools view -bS - > ./bam/sort_rmSX.bam
  # 或在call methylation步骤中，设置参数 “-x”（或“--rm-SX”）可移除被标记为 XS:i:1 的read
  rm merge.bam

  samtools index ./bam/sort_rmSX.bam

  samtools idxstats ./bam/sort_rmSX.bam |cut -f 1 |grep -e "[A-Za-z1-9]" |while read CHR
  do
    samtools view ./bam/sort_rmSX.bam $CHR -h | samtools view -Sb - > ${CHR}.bam
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/miniconda3/envs/py27/bin/python \
    /hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/BSseeker/bs_seeker2-call_methylation.py \
    -i ${CHR}.bam \
    -d ${index_position} \
    -o $CHR \
    --rm-CCGG \
    --rm-overlap
    rm ${CHR}.bam
    rm ${CHR}.bam_sorted.bam
    rm ${CHR}.bam_sorted.bam.bai
  done

  mv *.call_methylation_log ./log
  zcat  *.wig.gz > ./align_result/merge.wig
  rm *.wig.gz
  zcat *.ATCGmap.gz > ./align_result/merge.ATCGmap
  rm *.ATCGmap.gz
  zcat *.CGmap.gz > ./align_result/merge.CGmap
  rm *.CGmap.gz
  mkdir call_snp
  cgmaptools snv -i ./align_result/merge.ATCGmap -m bayes -v ./call_snp/bayes.vcf -o ./call_snp/bayes.snv --bayes-dynamicP
  workdir=`pwd`
  ATCGmap=$workdir/align_result/merge.ATCGmap
  CGmap=$workdir/align_result/merge.CGmap
  WIZ=$workdir/align_result/merge.wig
  vcf=$workdir/call_snp/bayes.vcf
  snv=$workdir/call_snp/bayes.snv
  echo -e "$sampleid\t$ATCGmap\t$CGmap\t$WIZ\t$vcf\t$snv" >> ../result_path.manifest
  echo "$sampleid have analysis done"

cd ../
done
echo "\nFinished !"
