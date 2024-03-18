---

layout: post

title: "It is new"

---

 My content

# Software used for RNA-seq(dRNA-seq and cRNA-seq)

## Convert to slow5 and convert back

```shell
#convert fast5 files to slow5 files using 8 I/O processes
slow5tools f2s fast5_dir -d blow5_dir  -p 8
#Merge all the slow5 files in to a single file using 8 threads
slow5tools merge blow5_dir -o file.blow5 -t8
#remove the temporary directory
rm -rf  blow5_dir
#Split the single SLOW5 file into multiple SLOW5 files such that each file has 4000 reads
slow5tools split file.blow5 -d blow5_dir -r 4000
#Now convert to FAST5 using using 8 I/O processes
slow5tools s2f blow5_dir -d fast5  -p 8
```

## Basecall

```shell
guppy_basecaller -c the_cfg_for_rna -i /fast5_files/ -s outputdir/
```

## Fastq data management

```shell
pyBioTools Fastq Filter 
-i theInputDocument
-o theOutputFile
#######some selectable options
--remove_duplicates 
--min_len 
--min_qual 
```

## Data quality assessment

**NanoPlot**: A tool uses fasta\q document from guppy. 

```shell
Nanoplot --summary summarydocumentbyguppy.txt --loglength -o output_dictionary
```

## Mapping and alignment

**gffread**

```shell
gffread genome.gtf -o- > genome.gff3
gffread -F -w transcriptome.fa -g genome.fa genome.gff3 
```

##### mapping:*: minimap2*

```shell
minimap2 -ax map-ont -splice -uf -k14 -t 4 -p 0 -N 10 the_reference_transcriptome.fa\
basecall_by_guppy.fastq > minimap_transcriptome.sam
```

**sort&index**

```shell
samtools sort -@ 4 -O bam -o xxx.bam xxx.sam 
#this step can automatically transform .sam  to /.bam
samtools index xxx.bam
samtools faidx xxx.fna
```

**show the quality of mapping**

```shell
samtools flagstat xxx.bam > xxx.txt
```

## Estimating Transcripts Count

**NanoCount**

```shell
NanoCount -i aligned_reads.bam -o transcript_counts.tsv
```

## Estimate poly(A) length(dRNA-seq)

**nanopolish-polya**

```shell
#### index
nanopolish index reads.fq --slow5 signals.blow5
#### poly(A) estimate
nanopolish polya --reads data_after_guppy.fastq --bam data_after_map.bam \
--genome reference.fas > ./polya_results.tsv
#### filter 'PASS'
grep 'PASS' polya_results.tsv > polya_results.pass_only.tsv
head -1 polya_results.tsv > header.tsv
cat header.tsv polya_results.pass_only.tsv > result.tsv
```

## RNA modification detection(dRNA-seq)

```shell
#### One can create a new conda environment for install nanocompore 
##### create environment with python 3.7
conda create -n nanocompore python=3.7
conda install -c bioconda nanocompore
#### data preparation
nanopolish index -s {sequencing_summary.txt} -d {raw_fast5_dir} {basecalled_fastq}

f5c_x86_64_linux_cuda eventalign --slow5 ../slow5/file.blow5 -b ../mapping/WT-1_transcript.bam -g ../../../../transcriptome.fa -r ../WT-1.fastq --rna > event_result.tsv

nanopolish eventalign --reads {basecalled_fastq} \
                      --bam {aligned_reads_bam}  \
                      --genome {transcriptome_fasta} \
                      --print-read-names \
                      --scale-events     \
                      --samples >{eventalign_reads_tsv}
nanocompore eventalign_collapse -t 6 \
                                -i {eventalign_reads_tsv} \
                                -o {eventalign_collapsed_reads_tsv}                      
#### start calculate
nanocompore sampcomp \
    --file_list1 ./data/S1_R1.tsv,./data/S1_R2.tsv \ #### R1/R2 is duplication which is recommanded but not necessary
    --file_list2 ./data/S2_R1.tsv,./data/S2_R2.tsv \
    --label1 S1 \
    --label2 S2 \
    --fasta ./reference/ref.fa \
    --outpath ./results
```

```shell
#### followed analysis needs python
from nanocompore.SampCompDB import SampCompDB, jhelp
db = SampCompDB (db_fn = "results/simulated_SampComp.db",
    fasta_fn = "references/simulated/ref.fa")
#### the path of db_fn and fasta_fn have better be the relative path.
#### Print general metadata information
print (db)
#### Prit list of references containing valid data
print (db.ref_id_list)


#### Generate text reports
#### It can create three types of files: save_report / save_shift_stats / save_to_bed
####### save_report : containing all the statistical results for all the positions
db.save_report (output_fn="./results/simulated_report.tsv")
#downstream edition
awk '{if($7!="nan"&&$7!="1.0")print}' simulated_report.tsv > simulated_report_final.tsv

####### save_shift_stats : Save the mean, median and sd intensity and dwell time for each condition and for each position.
db.save_shift_stats (output_fn="./results/simulated_shift.tsv")

####### save_to_bed
db.save_to_bed (output_fn="./results/simulated_sig_positions.bed")
```

## Alternative Splicing

**local AS events**

```shell
#Generation of local alternative splicing events
##This step aims to find all especial splicing events in reference genome by annotation file
 suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioe -e <list-of-events>({SE,SS,MX,RI,FL})
##<event-type>: correspond to the two letter code of the event from the following list.
##SE: Skipping Exon
##A5: Alternative 5' Splice Site
##A3: Alternative 3' Splice Site
##MX: Mutually Exclusive Exon
##RI: Retained Intron
##AF: Alternative First Exon
##AL: Alternative Last Exon
### If annotation file is not downloaded from Ensembl or Gencode (e.g.: RefSeq and UCSC genes), 
--pool-genes
# the input file must be gtf. IF you are gff, you need to convert it by gffread
#You can put all event files together
awk '
FNR==1 && NR!=1 { while (/^<header>/) getline; }
1 {print}
' *.ioe > allevents.ioe

#Combining multiple expression files
suppa.py joinFiles -f tpm -i sample1.tpm sample2.tpm sample3.tpm -o all_samples_tpms 

#PSI calculation for Transcripts and Events
suppa.py psiPerEvent --ioe-file <ioe-file> --expression-file <expression-file> -o <output-file> 
## the expression file's transcript name should be in accordance with annotation file
sed -i "s/^/rna-/g" virC_transcript_counts_TPM.tsv.tpm
#Differential splicing analysis
## You should have replication to do this step
## If not, this step can generate a file (not the required file for differential splicing analysis) which match event.
## You can analysis by deltaPSI(PSI(condition2)-PSI(condition1)) which can show the expression of different transcript.
suppa.py diffSplice --method <empirical> --input <ioe-file> --psi <Cond1.psi> <Cond2.psi> --tpm <Cond1_expression-file> <Cond2_expression-file> --area <1000> --lower-bound <0.05> -gc -o <output-file>
```

**transcript events**

```shell
#Generation of local alternative splicing events
##This step aims to find all especial splicing events in reference genome by annotation file
suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioi 

#PSI calculation for Transcripts and Events
suppa.py psiPerIsoform -g <gtf-file> -e <expression-file> -o <output-file>

#Differential splicing analysis
## You should have replication to do this step
suppa.py diffSplice --method <empirical> --input <ioi-file> --psi <Cond1.psi> <Cond2.psi> --tpm <Cond1_expression-file> <Cond2_expression-file> --area <1000> --lower-bound <0.05> -pa -gc -o <output-file>
```

**Cluster analysis**

```shell
#You must have replication in this step
suppa.py clusterEvents --dpsi <dpsi-file> --psivec <psivec-file> --sig-threshold <0.05> --eps <0.05> --min-pts <20> --groups <1-3,4-6> -o <output-file> 
```

