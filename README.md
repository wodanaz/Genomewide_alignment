# Genomewide alignment of high quality genomes (using sea urchins as example)
I will use Cactus to create a genome wide alignment


# First Step: Create an index/database and model repeat families of each genome with RepeatModeler


```bash

conda activate repeatmasker

####################### L. var

nano repeatmodel_lvar.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem 70000
BuildDatabase -name Lvar Lvar.fasta  # to build database
RepeatModeler -database Lvar -pa 23


# Run

sbatch repeatmodel_lvar.sh

####################### H. tub


nano repeatmodel_htub.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem 35000
BuildDatabase -name Htub Htub.fasta  # to build database
RepeatModeler -database Htub -pa 23


sbatch repeatmodel_htub.sh

####################### H. ery

nano repeatmodel_hery.sh 
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem 35000
BuildDatabase -name Hery Hery.fasta  # to build database
RepeatModeler -database Hery -pa 23



# run

sbatch repeatmodel_hery.sh 


#############################
# to check the memory allowance 
sacct -j ID --format=JobID,JobName,ReqMem,MaxRSS,Elapsed  # RAM requested/used!!


```

# Second Step: Soft-masking genomes with RepeatMasker using output from RepeatModeler.

```bash

####################### H. ery


mkdir Hery_mask_custom

nano repeatmasker_hery.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
RepeatMasker -engine ncbi -pa 23 -s -lib Hery-families.fa -gff -dir Hery_mask_custom -xsmall Hery.fasta

# run
sbatch repeatmasker_hery.sh


####################### H. tub

mkdir Htub_mask_custom

nano repeatmasker_htub.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
RepeatMasker -engine ncbi -pa 23 -s -lib Htub-families.fa -gff -dir Htub_mask_custom -xsmall Htub.fasta

# Run

sbatch repeatmasker_htub.sh

########################## L. var

mkdir Lvar_mask_custom

nano repeatmasker_lvar.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
RepeatMasker -engine ncbi -pa 23 -s -lib Lvar-families.fa -gff -dir Lvar_mask_custom -xsmall Lvar.fasta


# Run

sbatch repeatmasker_lvar.sh


```

# Third Step: Initialize Progressive Cactus to make genome wide alignments

Create urchin_seqfile.txt file with masked genomes.
The \* simbolizes the genomes are reference quality

```bash
nano urchin_seqfile.txt
  # Sequence data for progressive alignment of 3 genomes
  # Lv, He and Ht are flagged as good assemblies.
  # all will be used as an outgroup species.
(Lv:1,(He:0.2,Ht:0.2):0.8);
*He /data/wraycompute/alejo/PS_tests/Genome_alignments/masking_genomes/Hery.masked.fasta
*Ht /data/wraycompute/alejo/PS_tests/Genome_alignments/masking_genomes/Htub.masked.fasta
*Lv /data/wraycompute/alejo/PS_tests/Genome_alignments/masking_genomes/Lvar.masked.fasta


```




```bash
conda activate cactus
module load ucsc
mkdir urchins_wkdir

nano cactus_fastas2hal.sh
#!/usr/bin/env bash
#SBATCH -J maker_male
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 4
#SBATCH -c 24
#SBATCH --mem=40G
cactus --workDir urchins_wkdir/ --maxCores 96 --maxMemory 150G jobStore_urchin urchin_seqfile.txt urchins.hal --binariesMode local



sbatch cactus_fastas2hal.sh

```

Next, we need to convert HAL to MAF for each chromosome:


```bash
nano cactus_hal2maf.sh
#!/usr/bin/env bash
#SBATCH -J hal2maf
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH --mem=4G
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; do
hal2maf --noAncestors urchins.hal $chr.urchins.maf --refGenome He --refSequence $chr
done



sbatch cactus_hal2maf.sh

```


To extract features, we need to extract fasta using msa_split but first we need to proceess a bed file into feature data per chromosome

To create features

```bash
mkdir features
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; 
	do grep -w $chr Hery_adult_peaks.bed  | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done
```

Use msa_split to extract fasta files using the genome of reference.

```bash


awk '/^>chr/ {OUT=substr($0,2) ".fa";print " ">OUT}; OUT{print >OUT}' Hery_genome.fasta.masked

# Obtain a genome size file
module load samtools
samtools faidx Hery_genome.fasta.masked
awk '{print $1 "\t" $2 }' Hery_genome.fasta.masked.fai > sizes.genome




mkdir query


nano maf2fasta.sh
#!/usr/bin/env bash
#SBATCH -J maf2fasta
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY;
do msa_split $chr.urchins.maf --refseq $chr.fa --gap-strip ANY -q --in-format MAF --features features/$chr.feat.bed --for-features --out-root query/$chr; 
done



sbatch maf2fasta.sh
```









