# Genomewide_alignment
I will use Cactus to create a genome wide alignment


First, the genomes to align should be soft masked with RepeatMasker.


```bash
nano repeatmodel_lvar.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH --mem 30000
conda activate repeatmasker
BuildDatabase -name Lvar Lvar.fasta  # to build database
RepeatModeler -database Lvar -pa 23



sbatch repeatmodel_lvar.sh




# to check the memory allowance 
sacct -j ID --format=JobID,JobName,ReqMem,MaxRSS,Elapsed  # RAM requested/used!!


```

# Second Step: Masking genome with RepeatMasker using output from RepeatModeler.

```bash
nano repeatmasker_lvar.sh

#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
conda activate repeatmasker
./RepeatMasker -engine ncbi -pa 23 -s -lib ./Lvar_repeat_library.fa -gff -dir Lvar_mask_custom -xsmall /data/wraycompute/phil/urchin_genome/assemblies/scaffolds/Lvar_genome.fasta

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



