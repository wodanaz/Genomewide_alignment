# Genomewide_alignment
I will use Cactus to create a genome wide alignment




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



