# BSseq-WDL

# BSseq-WDl流程介绍

## 背景

BSseq-WDl 是我根据 BS-seeker2 + CGmapTools 为基础搭建的甲基化分析主要流程。说是主要，其实就是输出到甲基化结果文件（甲基化位点和甲基化率）。输出文件

> 可以类比为RNA-seq输出了表达矩阵，但是后续的差异，可视化这些我都没开始弄。原因主要是工作是做单细胞和空间组的。而且表观基因组学一直是我认为组学里面最难的。所以如果想做的完善，肯定要投入很多时间和精力，我只是业余时间搞的，后续有时间再搞。但是输出表达矩阵应该对很多分析是够用的。

## 该流程的优点

说这个流程的优点，其实就是说wdl的优点，相比于之前用shell写的。我觉得wdl最大优点就是：
- 彻底不用管文件的归档和存放，文件的输出和输入实现流程化控制
- 任务的并行运行，以及多个任务的多次调用比较方便
- 搭建虽然辛苦，但是维护起来比较简单

但是说完优点，也要说缺点，真的也是有点难以忍受：
- 没有积极活跃的社区，用户少，意味着问题很大概率要自己debug，让人崩溃
- workflow，call，task 三者都要一个input，再加上一个output。一直让人来回输入那几个变量（因为要变量承上启下，还要声明），繁琐且容易报错，让人没有感觉到一点代码的简洁之美。
- 文件放在特定的目录下，层级太多，很难找到，必须搭配输出的json文件
- if 语句中竟然还要评估否条件下的表达式，离谱！

其他的缺点暂时想不到了，如果一句话概括就是：

> 过程很痛苦，结果很美好，可能是流程控制语言的通病吧

## 流程空间流程

用wdl编写流程，需要cromwell.jar才能运行，输入文件有两个：

- 参数输入文件：*.json
- 样本名称和fq路径表格文件：*.tsv

输出文件需要添加 -m 参数，详情请见**运行**部分，输出的json文件：
包括了：

- 过滤后的read
- mapping后的用于call methylation的bam文件
- 甲基化结果文件`*.wig`, `*.ATCGmap`,`*.CGmap`
- call snp的结果文件`*.snv`,`*.vcf`

其中`*.ATCGmap`,`*.CGmap`是主要输出结果，前者代表是多种类型的甲基化位点结果（CHH，CHG，CG），后者只有CG类型的甲基化位点信息。结果是表格形式的，分别存储有位点，绝对甲基化reads数和相对的甲基化率等数值。

整体流程采用并行的方法，最大效率提高了运行速度和准确性：

![image-20220215091200415](https://pic-1259340288.cos.ap-guangzhou.myqcloud.com/img/202202150912001644887520Ovne4Bimage-20220215091200415.png)

## 环境要求

本流程运行需要以下环境：

- 需要python 2.7（pysam）
- 需要java 1.8.0 的环境
- 环境中要有samtools，fastp，[CGmapTools](https://cgmaptools.github.io)环境变量
- 需要下载[BSseeker2](https://github.com/BSSeeker/BSseeker2)

关于BSseeker和CGmapTools的安装示例问题，可以参考后文

## 流程运行

### 软件下载
下载Github仓库，解压，在联网环境下，运行`download_cromwell.sh`脚本，下载`cromwell-58.jar`,如果已经有`cromwell`其他版本的可以不用下载，可以自行尝试是否版本兼容。

```bash
bash download_cromwell.sh
```

### 配制文件

输入文件有两个，一个控制参数的`json`文件，一个是记录样本名称和fq序列文件路径的`tsv`表格.

### tsv 表格

制表符分割的文件，有三列，分别是样本名称，read1的路径，read2的路径。

Example:

```tsv
sample_A	./test/sample_A.test.1k.R1.fq.gz	./test/sample_A.test.1k.R2.fq.gz
sample_B	./test/sample_B.test.1k.R1.fq.gz	./test/sample_B.test.1k.R2.fq.gz
```

> 暂不考虑单端测序和生物学重复，前者是用不到，后者是还没做差异

### json文件

需要一些参数才能运行：

```json
{
    "bsseq.genome_index" : "error_path_just_for_test/test.fa_bowtie2/",
    "bsseq.python_file" : "wangjiaxuan/biosoft/miniconda3/envs/py27/bin/python",
    "bsseq.bsseeker_util_dir" : "wangjiaxuan/biosoft/BSseeker",
    "bsseq.refer_fa" : "./test/test.fa",
    "bsseq.alig_software" : "bowtie2",
    "bsseq.fastq_table" : "./test/input_sample_fq.tsv",
    "bsseq.step02_fq_qc" : "true",
    "bsseq.step03_bsseeker_align" : "true"
}
```
其中`genome_index`是基因组索引文件，这个是要BSseeker构建的。如果之前没有构建过，这个参数可以不写或者给了错误的地址。**流程都会自动根据参数`refer_fa`来进行基因组索引构建，只是需要的时间和内存都会比较大。

另外参数`refer_fa`是必须参数，除此以外，必须的参数还有python2.7路径的`python_file`，BSseeker的脚本目录`bsseeker_util_dir`，输入上文提到的tsv表格的`fastq_table`，以及比对软件`alig_software`。

### 运行

```bash
java -jar WDL/WDL/cromwell-58.jar run WDL/bsseq.wdl -i bsseq.input.json -m bsseq.output.json
# -m 参数是输出文件，可以不写，但是不写后续分析很麻烦
```

### 结果解读

结果可能有点疑惑，但这也是wdl的特性，输出文件都会在cromwell-executions文件里，里面根据task来进行分文件，又进一步细分为input和execution。但是看`bsseq.output.json`就不会疑惑了。

有用的信息，我都output了

Example：
```
  "outputs": {
    "bsseq.bsseeker_call_methylation_bam_bai": ["call-bsseeker_align_mergeBam/shard-0/execution/sort_rmSX.bam.bai", "call-bsseeker_align_mergeBam/shard-1/execution/sort_rmSX.bam.bai"],
    "bsseq.bsseeker_call_methylation_bam": ["call-bsseeker_align_mergeBam/shard-0/execution/sort_rmSX.bam", "call-bsseeker_align_mergeBam/shard-1/execution/sort_rmSX.bam"],
    "bsseq.sample": ["sample_A", "sample_B"],
    "bsseq.fastp_filter_report_json": ["call-bsseeker_fastq_qc/shard-0/execution/fastp.json", "call-bsseeker_fastq_qc/shard-1/execution/fastp.json"],
    "bsseq.fastp_filter_fq1": ["call-bsseeker_fastq_qc/shard-0/execution/sample_A.filter.R1.fq.gz", "call-bsseeker_fastq_qc/shard-1/execution/sample_B.filter.R1.fq.gz"],
    "bsseq.bsseeker_methy_CGmap": ["call-bsseeker_call_methy_merge/shard-0/execution/merge.CGmap", "call-bsseeker_call_methy_merge/shard-1/execution/merge.CGmap"],
    "bsseq.bsseeker_index": "call-bsseeker_build_index/execution/test.fa_bowtie2",
    "bsseq.bsseeker_methy_wig": ["call-bsseeker_call_methy_merge/shard-0/execution/merge.wig", "call-bsseeker_call_methy_merge/shard-1/execution/merge.wig"],
    "bsseq.cgmaptools_methy_snv": ["call-cgmaptools_call_snp/shard-0/execution/bayes.snv", "call-cgmaptools_call_snp/shard-1/execution/bayes.snv"],
    "bsseq.fastp_filter_fq2": ["call-bsseeker_fastq_qc/shard-0/execution/sample_A.filter.R2.fq.gz", "call-bsseeker_fastq_qc/shard-1/execution/sample_B.filter.R2.fq.gz"],
    "bsseq.cgmaptools_methy_vcf": ["call-cgmaptools_call_snp/shard-0/execution/bayes.vcf", "call-cgmaptools_call_snp/shard-1/execution/bayes.vcf"],
    "bsseq.bsseeker_methy_ATCGmap": ["call-bsseeker_call_methy_merge/shard-0/execution/merge.ATCGmap", "call-bsseeker_call_methy_merge/shard-1/execution/merge.ATCGmap"],
    "bsseq.fastp_filter_report_html": ["call-bsseeker_fastq_qc/shard-0/execution/fastp.html", "call-bsseeker_fastq_qc/shard-1/execution/fastp.html"]
  }
```

## 软件安装

### 安装BSseeker2
按照BSseeker的[网址](http://pellegrini-legacy.mcdb.ucla.edu/bs_seeker2/)上的介绍安装和使用
简单点话就是
```bash
git clone https://github.com/BSSeeker/BSseeker2
cd BSseeker2
```
就可以看到其中的py脚本，但是亲测python 3.8不能用，只能用conda创建一个python2.7的环境，同时需要环境中有bowtie, bowtie2, soap三个软件中任意一个都可以。python环境也要有`pysam`模块。

### 安装CGmapTools

```bash
git clone https://github.com/guoweilong/cgmaptools
bash install.sh
# add cgmaptools to your system PATH if necessary
# open ~/.bashrc
vi ~/.bashrc
# add the following lines to the end of ~/.bashrc
export PATH=/path/to/cgmaptools/:$PATH
# then, souce your ~/.bashrc
source ~/.bashrc
```

还需要在python2.7的环境中安装pysam和scipy
```bash
conda activate py27
pip install pysam
pip install scipy
```

## 资料来源

资料主要来源黄湘仪的教程——[18.10.04亚硫酸盐测序数据基本分析流程代码及软件](https://guoweilong.github.io/BSseeker2_CGmapTools_tutorial_181014.pdf)，当然也感谢BSseeker2 和CGmapTools开发者农大的[郭伟龙](http://guoweilong.github.io/)博士（中国农业大学）和朱平博士（天津血研所/北京大学）等人合作开发了DNA甲基化数据分析工具包：[CGmapTools](https://cgmaptools.github.io/)，包含40个独立的命令行工具，提供了全面的DNA甲基化数据分析和可视化功能。目前该工作已经在生物信息学著名期刊**Bioinformatics**上以original paper形式发表。

> Weilong Guo *#, Ping Zhu* , et al. (2017), [CGmapTools improves the precision of heterozygous SNV calls and supports allele-specific methylation detection and visualization in bisulfite-sequencing data](https://doi.org/10.1093/bioinformatics/btx595).

**再次感谢各位大佬创建这么好用的软件**