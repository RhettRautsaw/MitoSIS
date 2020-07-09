# BarcodeIDChecker

This script is a wrapper for many other programs and is designed extract the mitochondrial genome from raw or trimmed sequencing data and blast the resulting protein coding sequences or barcoding genes (e.g., CYTB, COX1, ND4, etc.) to check for sample mislabeling. 

**IN DEVELOPMENT**::The second part of this script will use the assembled mitogenomes to infer a phylogeny using HOMBLOCKS and IQTREE. 

## Input
**--assemble**

- **-f1, -f2**: Fastq files containing the raw/trimmed reads 

- **-r**: A fasta file containing reference sequences to map reads to prior to assembly and blast to following assembly. 
	- I recommend downloading **any** and **all** mitochondrial sequences/genomes for your clade from GenBank.
	- For example, since I study pitvipers I would download all mitochondrial sequences for snakes.  

**--phylo** (IN DEVELOPMENT)

- **-s**: A txt list containing the names of all samples on which assemble has been run. 

## Pipeline

**--assemble**

1. Map fastq reads to reference fasta using bwa 
2. Keep reads that successfully mapped using samtools 
3. Assemble mitogenome using MitoZ 
	- If fail, use best mitogenome and MITGARD to assemble mitogenome
4. Extract protein coding/barcoding genes 
5. Blast genes to reference fasta and determine top blast hit for taxa ID 

**--phylo**

1. Align mitogenomes using HomBlocks 
2. Infer phylogeny using IQTREE 

## Installation
**System Requirement**

- Linux
 
**Dependencies**

```
conda create -n mitgard_env \
biopython=1.69 \
blast=2.2.31 \
bowtie2=2.3.5 \
bwa=0.7.12 \
circos=0.69 \
ete3=3.0.0b35 \
gnu-parallel \
hmmer=3.1b2 \
infernal=1.1.1 \
iqtree \
libgd=2.2.4 \
minimap2=2.17 \
openjdk \
perl-bioperl \
perl-clone \
perl-list-moreutils \
perl-params-validate \
python=3.6.0 \
samtools=1.9 \
spades=3.13.1 \
tbl2asn \
trinity=2.8.5

git clone https://github.com/pedronachtigall/MITGARD.git
echo "export PATH=$PATH:path/to/MITGARD/bin/" >> ~/.bash_profile

git clone https://github.com/linzhi2013/MitoZ.git
tar -jxvf MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2
echo "export PATH=$PATH:path/to/release_MitoZ_v2.4-alpha/" >> ~/.bash_profile
```


## Example
```
BarcodeIDChecker_v1.py --assemble -f1 {}_F.fq -f2 {}_R.fq -r mito_ref.fasta -o {}

BarcodeIDChecker_v1.py --phylo -s samples.txt 
```


**Example PBS Script**

```
#PBS -N BarcodeIDChecker
#PBS -l select=20:ncpus=16:mem=100gb,walltime=72:00:00
#PBS -j oe
#PBS -M rrautsa@clemson.edu
#PBS -m abe

source .bash_profile

module load anaconda3
source activate mitgard_env

cd $PBS_O_WORKDIR

bwa index mito_ref.fasta
makeblastdb -in mito_ref.fasta -dbtype nucl

parallel -a samples.txt --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
	module load anaconda3
	source activate mitgard_env
	cd /path/to/workdir
	mkdir {}
	cd {}
	BarcodeIDChecker_v1.py --assemble -r ../mito_ref.fasta \
	-f1 /path/to/workdir/{}/02_trimmed/{}_F.fq.gz \
	-f2 /path/to/workdir/{}/02_trimmed/{}_R.fq.gz \
	-o {} -t 16"

BarcodeIDChecker_v1.py --phylo -s samples.txt
```

## Output
Inside the `BarcodeIDChecker_results` folder, you will find a tab-delimited file with the top 10 results for each successfully annotated mitochondrial gene for you to check your species ID including percent identity.

- `{}_blast.tab`

In addition, you will also find a `{}.result` folder that contains:

- The mitogenome `{}.fasta`
- The cds regions `{}.cds`
- A mitogenome plot `{}.circos.*`
- A genbank file `{}.mitoscaf.fa.gbf`

If `MitoZ` failed and `MITGARD` was required, you will also find the `best_reference.*` files which identify what reference sequence was used to create the mitogenome. The `{}_F.fq` and `{}_R.fq` files are the fastq reads which mapped to the provided reference sequences via `bwa`.

## Cite
- https://github.com/reptilerhett/BarcodeIDChecker
- MITGARD
- MitoZ
