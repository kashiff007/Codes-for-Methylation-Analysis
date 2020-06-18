## DNA methylation pipeline

![alt text](https://github.com/kashiff007/Codes-for-Methylation-Analysis/blob/master/DNA_methylation_pipeline.png)

Reads obtained from bisulfite sequencing in fastq format.

#### Filtering reads on the basis of quality
Reads are filteres to desired quality. I have used [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html) for this purpose. There are many commands in FASTX for different purpose, I am showing here an example to trim last 10 bases of each reads which has bad quality at the 3' end [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score#:~:text=A%20Phred%20quality%20score%20is,in%20the%20Human%20Genome%20Project) less than 25. Each read has a size of 150 base.

```{r, engine='bash', count_lines}
fastx_trimmer -i input.fastq -f 140 -o trimmed.fastq
```
#### Mapping BS-seq reads on genome through Bismark

Paired-end reads
```{r, engine='bash', count_lines}
~/Software/Bismark-0.22.3/bismark Genome_TAIR10/Original_Genome/ --multicore 4 --bowtie2 -p 10 -1 Methylation_reads/zr2939_16_R1.fq.gz -2 Methylation_reads/zr2939_16_R2.fq.gz --bam -o Col-0_SA187_Rep1
```
Single-end reads
```{r, engine='bash', count_lines}
~/Software/Bismark-0.22.3/bismark --multicore 2  -p 4 --bam -o Symb_bismark_Rep1 Genome_TAIR10/Original_Genome/ ../Maha_ChiP_RNA/S1_1.fq
```

Remove Duplication
