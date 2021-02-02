# Genomewide_alignment
I will use Cactus to create a genome wide alignment


First, the genomes to align should be soft masked with RepeatMasker.


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
#conda activate repeatmasker
BuildDatabase -name Htub Htub.fasta  # to build database
RepeatModeler -database Htub -pa 23


sbatch repeatmodel_htub.sh

####################### H. ery

nano repeatmodel_hery.sh 
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem 35000
conda activate repeatmasker
BuildDatabase -name Hery Hery.fasta  # to build database
RepeatModeler -database Hery -pa 23



# run

sbatch repeatmodel_hery.sh 


#############################
# to check the memory allowance 
sacct -j ID --format=JobID,JobName,ReqMem,MaxRSS,Elapsed  # RAM requested/used!!


```

# Second Step: Masking genome with RepeatMasker using output from RepeatModeler.

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

Create urchin_seqfile.txt file with masked genomes.
the \* simbolizes the genomes are reference quality

```bash
nano urchin_seqfile.txt
  # Sequence data for progressive alignment of 4 genomes
  # human, chimp and gorilla are flagged as good assemblies.
  # since orang isn't, it will not be used as an outgroup species.
 (((human:0.006,chimp:0.006667):0.0022,gorilla:0.008825):0.0096,orang:0.01831);
 *human /data/genomes/human/human.fa
 *chimp /data/genomes/chimp/
 *gorilla /data/genomes/gorilla/gorilla.fa
 orang /cluster/home/data/orang/


```




```bash
conda activate cactus
module load ucsc

nano cactus.sh
#!/usr/bin/env bash
#SBATCH -J maker_male
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 4
#SBATCH -c 24
#SBATCH --mem=40G

cactus --workDir /data/wraycompute/alejo/PS_tests/Genome_alignments/urchins/urchin2wkdir --maxCores 96 --maxMemory 150G jobStore_urchin2 urchin_seqfile.txt urchins2.hal --binariesMode local


```

To extract fasta files from a feature bed file:

1. we need to convert HAL to MAF for each chromosome


```bash
nano cactus_hal2maf.sh
#!/usr/bin/env bash
#SBATCH -J hal2maf
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -N 1
#SBATCH --mem=4G
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; do
hal2maf --noAncestors urchins2.hal $chr.urchins.maf --refGenome He --refSequence $chr ;
sed -ie 2d $chr.urchins.maf
done


sbatch cactus_hal2maf.sh

```


2. we need to extract fasta using msa_split

First, create features
```bash
mkdir features
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ; 
	do grep -w $chr Hery_adult_peaks.bed  | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done
```

Second, use msa_split to extraxt fasta files

```bash


awk '/^>chr/ {OUT=substr($0,2) ".fa";print " ">OUT}; OUT{print >OUT}' Hery_genome.fasta.masked


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









