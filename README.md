# **Project pipeline**

## **BAM pre-processing**

### **Sorting and indexing**
Sorting and indexing with samtool both for *Control.bam* and *Tumor.bam*

```bash
samtools sort input_files/Tumor.bam -o sorted_indexed_BAM/tumor.sorted.bam 

samtools index sorted_indexed_BAM/tumor.sorted.bam
```

flagstats for *control.sorted.bam :*

```
19720171 + 0 in total (QC-passed reads + QC-failed reads)
19708438 + 0 primary
11733 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
19669924 + 0 mapped (99.75% : N/A)
19658191 + 0 primary mapped (99.75% : N/A)
19708438 + 0 paired in sequencing
9854219 + 0 read1
9854219 + 0 read2
19576046 + 0 properly paired (99.33% : N/A)
19613806 + 0 with itself and mate mapped
44385 + 0 singletons (0.23% : N/A)
16828 + 0 with mate mapped to a different chr
10030 + 0 with mate mapped to a different chr (mapQ>=5)
```

flagstats for *tumor.sorted.bam :*

```
15039503 + 0 in total (QC-passed reads + QC-failed reads)
15029250 + 0 primary
10253 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
15033690 + 0 mapped (99.96% : N/A)
15023437 + 0 primary mapped (99.96% : N/A)
15029250 + 0 paired in sequencing
7514625 + 0 read1
7514625 + 0 read2
14979936 + 0 properly paired (99.67% : N/A)
15019614 + 0 with itself and mate mapped
3823 + 0 singletons (0.03% : N/A)
14070 + 0 with mate mapped to a different chr
7572 + 0 with mate mapped to a different chr (mapQ>=5)
```

Full statistics are in *tumor.stats.txt* and *control.stats.txt*

### **Filtering**

Both *control.sorted.bam* and *tumor.sorted.bam* have been filtered with samtools to meet the following criteria :

+ MAPQ over 30
+ Pass filters such as platform/vendor quality controls (flags)
+ being properly aligned according to the aligner

```bash
samtools view -q 30 -F 512 -f 2 -b control.sorted.bam > ../filtered_sorted_indexed_bam/filtered.control.sorted.bam


```

Doing so filtered out 23,11% of control reads (from 19.720.171 to 15.161.625) and 22.49% of tumor reads (from 15.039.503 to 11.656.539)

Filtered bam files were re-indexed.

### **Realignment around indels**

Realignment around indels was done with GATK using *human g1k v37* as reference genome and a list of high confidence SNPs in *.vcf* format from [the Broad Institute public data repository](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0).

First the two *.intervarls* files are generated :

```bash
java -jar ../../tools/GenomeAnalysisTK.jar -T RealignerTargetCreator
    -R ../annotations/human_g1k_v37.fasta
    -I ../filtered_sorted_indexed_bam/filtered.control.sorted.bam 
    -known ../annotations/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf
    -o control.realigner.intervals

java -jar ../../tools/GenomeAnalysisTK.jar -T RealignerTargetCreator
    -R ../annotations/human_g1k_v37.fasta
    -I ../filtered_sorted_indexed_bam/filtered.tumor.sorted.bam 
    -known ../annotations/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf
    -o tumor.realigner.intervals
```

once the two *.intervals* are created, the actual realignment can be performed :

```bash
java -jar ../../tools/GenomeAnalysisTK.jar -T IndelRealigner
-R ../annotations/human_g1k_v37.fasta
-I ../filtered_sorted_indexed_bam/filtered.control.sorted.bam
-known ../annotations/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf
-targetIntervals control.realigner.intervals
-o realigned_filtered_control_sorted.bam

java -jar ../../tools/GenomeAnalysisTK.jar -T IndelRealigner
-R ../annotations/human_g1k_v37.fasta
-I ../filtered_sorted_indexed_bam/filtered.tumor.sorted.bam
-known ../annotations/hg19_v0_1000G_phase1.snps.high_confidence.b37.vcf
-targetIntervals tumor.realigner.intervals
-o realigned_filtered_tumor_sorted.bam
```

GATK realigned 8305 reads (~5%) in the control bam and 6651 (~6%) in the tumor bam.

### **Duplicate removal**

Duplicate removal was conducted with samtools.

As first step, the bam files are sorted by read name (needed by samtools fixmate)

```bash
samtools sort -n ../realigned_filtered_sorted_indexed_bam/realigned_filtered_control_sorted.bam -o nsorted_realigned_filtered_control_sorted.bam

samtools sort -n ../realigned_filtered_sorted_indexed_bam/realigned_filtered_tumor_sorted.bam -o nsorted_realigned_filtered_tumor_sorted.bam
```

Add mate tags when not present

```bash
samtools fixmate -m nsorted_realigned_filtered_control_sorted.bam fixed_realigned_filtered_control_sorted.bam

samtools fixmate -m nsorted_realigned_filtered_tumor_sorted.bam fixed_realigned_filtered_tumor_sorted.bam
```

Now we can re-sort the bam by coordinates and run the actual duplicate removal,

re-sorting :

```bash
samtools sort fixed_realigned_filtered_control_sorted.bam > fixed_realigned_filtered_control_re-sorted.bam

samtools sort fixed_realigned_filtered_tumor_sorted.bam > fixed_realigned_filtered_tumor_re-sorted.bam
```

duplicate removal :

```
samtools markdup -sr fixed_realigned_filtered_control_re-sorted.bam dedup_realigned_filtered_control.bam

samtools markdup -sr fixed_realigned_filtered_tumor_re-sorted.bam dedup_realigned_filtered_tumor.bam
```

This removed 2.074.733 duplicate reads in the control bam (~14%) leading to a total of 13.086.892 control, and 1.379.650 duplicate reads in the tumor bam (~12%) leading to a total of 10.276.889 tumor reads.

To compare the results, duplicate removal was also done using picard :

```
java -jar ../../tools/picard.jar MarkDuplicates I= ../realigned_filtered_sorted_indexed_bam/realigned_filtered_control_sorted.bam O=dedup_picard_realigned_filtered_control.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=control_picard.log ASSUME_SORTED=true

java -jar ../../tools/picard.jar MarkDuplicates I= ../realigned_filtered_sorted_indexed_bam/realigned_filtered_tumor_sorted.bam O=dedup_picard_realigned_filtered_tumor.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=tumor_picard.log ASSUME_SORTED=true
```
Similarly to samtools, Picard removed 2.070.030 (~13%) duplicate reads from the control bam and 1.376.758 (~12%) from the tumor bam

### **Base quality score recalibration**