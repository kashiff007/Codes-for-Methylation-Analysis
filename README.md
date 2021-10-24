## DNA methylation pipeline

![alt text](https://github.com/kashiff007/Codes-for-Methylation-Analysis/blob/master/DNA_methylation_pipeline.png)

Reads obtained from bisulfite sequencing in fastq format.

#### Filtering reads on the basis of quality
Reads are filteres to desired quality. I have used [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html) for this purpose. There are many commands in FASTX for different purpose, I am showing here an example to trim last 10 bases of each reads which has bad quality at the 3' end [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score#:~:text=A%20Phred%20quality%20score%20is,in%20the%20Human%20Genome%20Project) less than 25. Each read has a size of 150 base.

```{r, engine='bash', count_lines}
fastx_trimmer -i input.fastq -f 140 -o trimmed.fastq
```
#### Mapping BS-seq reads on genome through Bismark

*Paired-end reads*
```{r, engine='bash', count_lines}
~/Software/Bismark-0.22.3/bismark Genome_TAIR10/Original_Genome/ --multicore 4 --bowtie2 -p 10 -1 condition_BS-seq_reads/condition1_R1.fq.gz -2 condition_BS-seq_reads/condition1_R2.fq.gz --bam -o Col-0_SA187_Rep1
```
*Single-end reads*
```{r, engine='bash', count_lines}
~/Software/Bismark-0.22.3/bismark --multicore 2  -p 4 --bam -o condition1_bismark_Rep1 Genome/Original_Genome/ ../condition_BS-seq_reads/condition1.fq
```

**Note: Do NOT sort the bam/sam file** 

#### Remove Duplicate reads from BAM file

Usually bismark map reads on basis of unique position, but still we have to remove the duplicates reads from PCR to avoid redundancy.
```{r, engine='bash', count_lines}
~/Software/Bismark-0.22.3/deduplicate_bismark --bam -p condition1_bismark_bt2_pe.bam
```

#### Calling DNA methylation from bismark methylation extractor

Now final bam file has all the information of mapped reads, we can use bam file as an input to obtain the DNA methylation for every location of Cs in the genome. 
```{r, engine='bash', count_lines}
~/Software/Bismark-0.22.3/bismark_methylation_extractor -p -o methylation_extractor_output --multicore 4 --comprehensive --merge_non_CpG --bedGraph --counts --scaffolds --CX --cytosine_report --CX --genome_folder .FULLPATH/Genome/Original_Genome/ condition1_bismark_bt2_pe_deduplicate.bam
```
**Note: Use full path for --genome_folder otherwise it will not generate genome-wide cytosine files**

This command generates 3 context files CpG/CG, CHG and CHH which has the information of each cytosine location for corresponding contexts. 
- CpG_context_file_name.text
- CHG_context_file_name.text
- CHH_context_file_name.text

The methylation extractor output looks like this (tab separated):
1. seq-ID
2. methylation state
3. chromosome
4. start position (= end position)
5. methylation call

Methylated cytosines will receive a '+' orientation, unmethylated cytosines will receive a '-'
orientation. 

       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ~~~   X   for methylated C in CHG context                      ~~~
       ~~~   x   for not methylated C CHG                             ~~~
       ~~~   H   for methylated C in CHH context                      ~~~
       ~~~   h   for not methylated C in CHH context                  ~~~
       ~~~   Z   for methylated C in CpG context                      ~~~
       ~~~   z   for not methylated C in CpG context                  ~~~
       ~~~   U   for methylated C in Unknown context (CN or CHN       ~~~
       ~~~   u   for not methylated C in Unknown context (CN or CHN)  ~~~
       ~~~   .   for any bases not involving cytosines                ~~~
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Examples for cytosines in CpG context:
```
HWUSI-EAS611_0006:3:1:1058:15806#0/1 - 6 91793279 z
HWUSI-EAS611_0006:3:1:1058:17564#0/1 + 8 122855484 Z
```
Examples for cytosines in CHG context:
```
HWUSI-EAS611_0006:3:1:1054:1405#0/1 - 7 89920171 x
HWUSI-EAS611_0006:3:1:1054:1405#0/1 + 7 89920172 X
```
Examples for cytosines in CHH context:
```
HWUSI-EAS611_0006:3:1:1054:1405#0/1 - 7 89920184 h 
```

Other files generated are:
* *.bismark.cov.gz
* *.bedGraph.gz
* *.CX_report.txt

Bismark.cov.gz has the information of coverage percentage and methylated and non-methylated read numberes for all cytosine.
```
Chr2    1006    1006    1.72413793103448        1       57
Chr2    1007    1007    2.04081632653061        1       48
Chr2    1009    1009    2.85714285714286        2       68
Chr2    1010    1010    3.33333333333333        2       58
Chr2    1012    1012    1.23456790123457        1       80
```
Bedgraph file has the information only about the coverge.
```
Chr2    1005    1006    1.72413793103448
Chr2    1006    1007    2.04081632653061
Chr2    1008    1009    2.85714285714286
Chr2    1009    1010    3.33333333333333
Chr2    1011    1012    1.23456790123457
```
Most important output file from bismark_methylation_extractor command is **CX_report.txt**. Because thsi file has the information of context, strand orientation and number of methylated and un-methylated reads.
```
Chr2    1001    -       0       0       CHH     CNN
Chr2    1006    +       1       57      CG      CGT
Chr2    1007    -       1       48      CG      CGA
Chr2    1009    +       2       68      CG      CGA
Chr2    1010    -       2       58      CG      CGA
Chr2    1012    +       1       80      CHH     CCA
Chr2    1013    +       4       88      CHG     CAG
Chr2    1015    -       2       83      CHG     CTG
```

Use this **CX_report.txt** file as an input for methylKit pipeline. But first convert the CX file in the format which methylKit accepts
```
awk '$6=="CG"{if($3=="+"&&($4>0&&$5>0 || $4>0&&$5==0))print$1"."$2"\t"$1"\t"$2"\tF\t"$4+$5"\t"$4/($4+$5)*100"\t"$5/($4+$5)*100; else if ($3=="-"&&($4>0&&$5>0 ||  $4>0&&$5==0)) print$1"."$2"\t"$1"\t"$2"\tR\t"$4+$5"\t"$4/($4+$5)*100"\t"$5/($4+$5)*100; else if ($3=="+"&&($4==0&&$5>0 || $4==0&&$5==0)) print$1"."$2"\t"$1"\t"$2"\tF\t"$4+$5"\t0\t0"; else if ($3=="-"&&($4==0&&$5>0 || $4==0&&$5==0)) print$1"."$2"\t"$1"\t"$2"\tR\t"$4+$5"\t0\t0"}' **.deduplicated.CX_report.txt > **.methylKit_input.txt 
```

The file for methylKit will look like:
```
##         chrBase   chr    base strand coverage freqC  freqT
## 1 chr21.9764539 chr21 9764539      R       12 25.00  75.00
## 2 chr21.9764513 chr21 9764513      R       12  0.00 100.00
## 3 chr21.9820622 chr21 9820622      F       13  0.00 100.00
## 4 chr21.9837545 chr21 9837545      F       11  0.00 100.00
## 5 chr21.9849022 chr21 9849022      F      124 72.58  27.42
```

Now we have files for methylation at per base level, we can use following code in R MethyKit.

```
library(methylKit)
library(tidyverse)
## Load methylKit files into list
file.list = list("H2-25-1_CG_methylKit_input.txt", "H2-25-2_CG_methylKit_input.txt", "H2-25-3_CG_methylKit_input.txt", "H2-25-5_CG_methylKit_input.txt", "H2-25-6_CG_methylKit_input.txt", "H2-25-7_CG_methylKit_input.txt", "H2-32-1_CG_methylKit_input.txt", "H2-32-2_CG_methylKit_input.txt", "H2-32-3_CG_methylKit_input.txt", "H2-32-5_CG_methylKit_input.txt", "H2-32-6_CG_methylKit_input.txt", "H2-32-7_CG_methylKit_input.txt")

## Read the files
H2_myobj=methRead(file.list, sample.id=list("H2_25_1","H2_25_2","H2_25_3","H2_25_5","H2_25_6","H2_25_7","H2_32_1","H2_32_2","H2_32_3","H2_32_5","H2_32_6", "H2_32_7" ), assembly="H2", treatment=c(0,0,0,0,0,0,1,1,1,1,1,1),header=FALSE, context="CpG", dbtype = NA)

## Filter each file based in read counts
H2_filtered.myobj=filterByCoverage(H2_myobj,lo.count=5,lo.perc=NULL,hi.count=NULL)

## Destand files: CG has two adjacent loacation, chnage it into one by adding read counts from C and G.
H2_filtered_myobj_meth = methylKit::unite(H2_filtered.myobj, destrand=TRUE, min.per.group=1L)

## change H2_filtered_myobj_meth into data frame
H2_filtered_myobj_meth_df <- data.frame(H2_filtered_myobj_meth)

## Replace all NA with 0
H2_filtered_myobj_meth_df[is.na(H2_filtered_myobj_meth_df)] <- 0

```

