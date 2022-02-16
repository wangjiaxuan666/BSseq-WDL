version 1.0

workflow bsseq {
    input {
        Boolean step02_fq_qc = true
        Boolean step03_bsseeker_align = true
        File python_file
        File bsseeker_util_dir
        String? genome_index
        File refer_fa
        String alig_software
        File fastq_table
        Array[Array[String]] sample_tsv = read_tsv(fastq_table)
        Int liner_num = 4
    }
    call bsseeker_build_index {
            input :
                python_file = python_file,
                run_index = bsseeker_util_dir + "/bs_seeker2-build.py",
                refer_fa = refer_fa,
                alig_software = alig_software,
                genome_index = genome_index,
                index_exist = if defined(genome_index) then "true" else "false"
    }

    scatter ( sample_sheet in sample_tsv ) {
        if (step02_fq_qc) {
            call bsseeker_fastq_qc {
                input :
                    sample_sheet = sample_sheet
            }
        }
    
        if ( step03_bsseeker_align ){
            call bsseeker_align_splitFq {
                input : 
                    sample_name = bsseeker_fastq_qc.sample_name,
                    filter_fq1 = bsseeker_fastq_qc.filter_fq1,
                    filter_fq2 = bsseeker_fastq_qc.filter_fq2
            }

            scatter (i in range(length(bsseeker_align_splitFq.split_read1))){
                call bsseeker_align_bam {
                    input :
                        python_file = python_file,
                        run_align = bsseeker_util_dir + "/bs_seeker2-align.py",
                        run_antisense_read = bsseeker_util_dir + "/Antisense.py",
                        refer_fa = refer_fa,
                        alig_software = alig_software,
                        liner_num = liner_num,
                        split_read1 = bsseeker_align_splitFq.split_read1[i],
                        split_read2 = bsseeker_align_splitFq.split_read2[i],
                        genome_db =  bsseeker_build_index.index
                        # genome_db = select_first([bsseeker_build_index.index,genome_index])
                }
            }
            call bsseeker_align_mergeBam {
                input:
                sample_name = bsseeker_fastq_qc.sample_name,
                bams = bsseeker_align_bam.bam
            }

            scatter (CHR in bsseeker_align_mergeBam.id) {
                call bsseeker_call_methy {
                    input : 
                    CHR = CHR,
                    python_file = python_file,
                    run_call_methy = bsseeker_util_dir + "/bs_seeker2-call_methylation.py",
                    genome_db =  bsseeker_build_index.index,
                    bs_bam = bsseeker_align_mergeBam.bs_bam
                }
            }

            call bsseeker_call_methy_merge {
                input :
                    sample_name = bsseeker_fastq_qc.sample_name,
                    wig = bsseeker_call_methy.wig,
                    ATCGmap = bsseeker_call_methy.ATCGmap,
                    CGmap = bsseeker_call_methy.CGmap
            }

            call cgmaptools_call_snp {
                input:
                sample_name = bsseeker_fastq_qc.sample_name,
                ATCGmap = bsseeker_call_methy_merge.ATCGmap
            }
        }
    }

    output {
        String? bsseeker_index  = bsseeker_build_index.index
        #Array[String?]? sample = select_first([bsseeker_fastq_qc.sample_name])
        Array[File?]? fastp_filter_fq1 = bsseeker_fastq_qc.filter_fq1
        Array[File?]? fastp_filter_fq2 = bsseeker_fastq_qc.filter_fq2
        #Array[File?]? fastp_filter_report_json = bsseeker_fastq_qc.filter_report_json
        #Array[File?]? fastp_filter_report_html = bsseeker_fastq_qc.filter_report_html
        Array[File?]? bsseeker_call_methylation_bam = bsseeker_align_mergeBam.bs_bam
        Array[File?]? bsseeker_call_methylation_bam_bai = bsseeker_align_mergeBam.bs_bai
        Array[File?]? bsseeker_methy_wig = bsseeker_call_methy_merge.wig
        Array[File?]? bsseeker_methy_ATCGmap = bsseeker_call_methy_merge.ATCGmap
        Array[File?]? bsseeker_methy_CGmap = bsseeker_call_methy_merge.CGmap
        Array[File?]? cgmaptools_methy_snv = cgmaptools_call_snp.snv
        Array[File?]? cgmaptools_methy_vcf = cgmaptools_call_snp.vcf
    }
}

task bsseeker_build_index {
    input {
        String python_file
        String run_index
        File refer_fa
        String refer_index_name = basename(refer_fa)
        String alig_software
        String index_exist
        String?  genome_index
    }
    command {
        if [ ${index_exist} == "false" ] || [ ! -e ${genome_index} ];then
        ${python_file} \
        ${run_index} \
        -f ${refer_fa} \
        --aligner=${alig_software} \
        -d ./
        else
        ln -s ${genome_index} ./
        fi
    }
    output {
        File? index = "${refer_index_name}" + "_" + "${alig_software}"
    }
}

task bsseeker_fastq_qc {
    input {
        Array[String] sample_sheet
        String sample_name = sample_sheet[0]
        String fq1 = sample_sheet[1]
        String fq2 = sample_sheet[2]
    }

    command {
        fastp \
        -i ${fq1} \
        -I ${fq2} \
        -o ${sample_name}.filter.R1.fq.gz \
        -O ${sample_name}.filter.R2.fq.gz
    }

    output {
        String sample_name = sample_name
        File filter_fq1 = "${sample_name}.filter.R1.fq.gz"
	    File filter_fq2 = "${sample_name}.filter.R2.fq.gz"
        File filter_report_json = "fastp.json"
        File filter_report_html = "fastp.html"
    }
}

task bsseeker_align_splitFq {
    input {
        String? sample_name
        File? filter_fq1 
        File? filter_fq2
    }
    command {
        zcat ${filter_fq1} | split -l 20000000 - R1_
        zcat ${filter_fq2} | split -l 20000000 - R2_
    }
    # 20,000,000 reads
    output {
        Array[String] split_read1 = glob("R1_*")
        Array[String] split_read2 = glob("R2_*")
    }
}

task bsseeker_align_bam {
    input {
        String python_file
        String run_align
        String run_antisense_read
        File refer_fa
        String? genome_db
        String alig_software 
        File split_read1
        File split_read2
        Int liner_num
        String? tmp_name1 = basename(split_read1)
        String? tmp_name2 = basename(split_read2)              
    }
    command {
        ${python_file} \
        ${run_align} \
        -i ${split_read1} \
        -g ${refer_fa} \
        -d ${genome_db}/../ \
        --aligner=${alig_software} \
        --bt2-p ${liner_num} \
        --bt2--end-to-end --bt2--very-sensitive --bt2--dovetail \
        --temp_dir=`pwd` \
        -m 0.1 \
        --XSteve \
        -o ${tmp_name1}.bam

        ${python_file} \
        ${run_antisense_read} \
        -i ${split_read2} \
        -o ${split_read2}.anti.fq
        
        ${python_file} \
        ${run_align} \
        -i ${split_read2}.anti.fq \
        -g ${refer_fa} \
        -d ${genome_db}/../ \
        --aligner=${alig_software} \
        --bt2-p ${liner_num} \
        --bt2--end-to-end --bt2--very-sensitive --bt2--dovetail \
        --temp_dir=`pwd` \
        -m 0.1 \
        --XSteve \
        -o ${tmp_name2}_anti.bam
    }
    output {
        Array[String] bam = glob("R[1-2]_[a-z][a-z]*.bam")
        Array[String] bam_log = glob("R[1-2]_[a-z][a-z]*.bam*log")
    }
}

task bsseeker_align_mergeBam {
    input {
        String? sample_name
        Array[Array[File]] bams
        Array[String] bams_l = flatten(bams)
    }

    command {
        samtools merge merge.bam ${sep=' ' bams_l}
        samtools sort merge.bam -o ${sample_name}_sort.bam
        samtools view -h ${sample_name}_sort.bam | grep -v 'XS:i:1' | samtools view -bS - > ${sample_name}_sort_rmSX.bam
        samtools index ${sample_name}_sort_rmSX.bam
        rm merge.bam 
        samtools idxstats ${sample_name}_sort_rmSX.bam |cut -f 1 |grep -e "[A-Za-z1-9]" > chr.id
    }
    output {
        File merge_sort_bam = "${sample_name}_sort.bam"
        File bs_bam = "${sample_name}_sort_rmSX.bam"
        File bs_bai = "${sample_name}_sort_rmSX.bam.bai"
        Array[String] id = read_lines("chr.id")
    }
}

task bsseeker_call_methy {
    input {
        String CHR
        String python_file
        String run_call_methy
        String? genome_db
        File bs_bam
    }

    command {
    samtools view ${bs_bam} ${CHR} -h | samtools view -Sb - > ${CHR}.bam
    ${python_file} \
    ${run_call_methy} \
    -i ${CHR}.bam \
    -d ${genome_db} \
    -o ${CHR} \
    --rm-CCGG \
    --rm-overlap
    rm ${CHR}.bam
    rm ${CHR}.bam_sorted.bam
    rm ${CHR}.bam_sorted.bam.bai
    }
    
    output {
        #File log = glob("*.log")
        File wig = "${CHR}.wig.gz"
        File ATCGmap = "${CHR}.ATCGmap.gz"
        File CGmap = "${CHR}.CGmap.gz"
    }
}

task bsseeker_call_methy_merge {
    input{
        String? sample_name
        Array[String] wig
        Array[String] ATCGmap
        Array[String] CGmap
    }

    command {
        zcat ${sep=' ' wig} > ${sample_name}.wig
        zcat ${sep=' ' ATCGmap} > ${sample_name}.ATCGmap
        zcat ${sep=' ' CGmap} > ${sample_name}.CGmap
    }

    output {
        File wig = "${sample_name}.wig"
        File ATCGmap = "${sample_name}.ATCGmap"
        File CGmap = "${sample_name}.CGmap"
    }
}

task cgmaptools_call_snp {
    input {
        String? sample_name
        File ATCGmap
    }
    command {
        cgmaptools snv \
        -i ${ATCGmap} \
        -m bayes \
        -v ${sample_name}_bayes.vcf \
        -o ${sample_name}_bayes.snv \
        --bayes-dynamicP
    }
    
    output {
        File snv = "${sample_name}_bayes.snv"
        File vcf = "${sample_name}_bayes.vcf"
    }
}