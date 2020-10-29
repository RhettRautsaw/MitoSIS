# BarcodeIDChecker
## Rhett M. Rautsaw

This script is a wrapper for [MitoZ](https://github.com/linzhi2013/MitoZ) and [MITGARD](https://github.com/pedronachtigall/MITGARD) and is designed extract vertebrate mitochondrial genomes from raw or trimmed sequencing data. It will then blast the resulting mitogenome or barcoding genes (e.g., CYTB, COX1, ND4, 16S, etc.) to check for sample mislabeling and produce a phylogeny using [IQ-TREE](http://www.iqtree.org/).

<br>

# Arguments

|       flag      |   description   |
|-----------------|-----------------|
| -h, --help      | Show this help message and exit. | 
| -f1, --fastq1   | fastq read pair 1 (forward). **No default setting.** |
| -f2, --fastq2   | fastq read pair 2 (reverse). **No default setting.** |
| -r, --reference | genbank database. **No default setting.** <br> *Recommend downloading all mitochondrial data for your clade of interest <br> e.g., snakes; [Genbank Example](https://www.ncbi.nlm.nih.gov/nuccore/?term=snakes%5Bporgn%5D+AND+mitochondrion%5Bfilter%5D) <br> Send to > Complete Record > Genbank* |
| -o, --output    | Prefix for output files. **Default is 'ZZZ'** |
| -c, --cpu       | Number of threads to be used in each step. **Default is 8** |
| -M, --memory    | Max memory for Trinity (see [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity) for format). **Default is '30G'** |
| --version       | Show program's version number and exit |
|<img width=200/> |<img width=500/>|

# Pipeline

1. Map fastq reads to reference fasta using bwa
2. Keep reads that successfully mapped using samtools
3. Assemble mitogenome using MitoZ
4. If MitoZ fails, identify the best reference sequence
5. Use MITGARD to assemble mitogenome
	- Use MitoZ to annotate mitogenome (if not already done)
6. Extract protein coding/barcoding genes
7. Blast mitogenome or genes to reference database
8. Export results and sequences
9. Align sequences and build phylogeny

![](BarcodeIDChecker_Flowchart.png)


# Installation
**System Requirement**
- Linux

**Conda Installation**
```
# Clone this Repository
git clone https://github.com/reptilerhett/BarcodeIDChecker.git
cd BarcodeIDChecker
echo "export PATH=\$PATH:$PWD" >> ~/.bash_profile

# Clone MITGARD Repository 
git clone https://github.com/pedronachtigall/MITGARD.git
echo "export PATH=\$PATH:$PWD/MITGARD/bin" >> ~/.bash_profile

# Clone MitoZ Repository
git clone https://github.com/linzhi2013/MitoZ.git
tar -jxvf MitoZ/version_2.4-alpha/release_MitoZ_v2.4-alpha.tar.bz2
echo "export PATH=\$PATH:$PWD/release_MitoZ_v2.4-alpha" >> ~/.bash_profile

# Make sure everythig has proper permissions and source your bash_profile
chmod -R 755 *
source ~/.bash_profile

# Create Conda Environment
conda env create -f barcode_env.yml
conda activate barcode_env

# Install dfply
pip install dfply

# Install Taxonomy Database for MitoZ
python
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
quit()

# YOU'RE READY TO GO
# Check if BarcodeIDChecker.py is in your path
BarcodeIDChecker.py -h
```


# Example
We recommend trimming your data first prior to running this program. Example trimming using [Trim-Galore](https://github.com/FelixKrueger/TrimGalore) shown below. Depending on whether you are working with DNA or RNA-Seq data, you may want to change the length/quality parameters.
```
# Trimming
trim_galore --paired --phred33 --length 30 -q 20 -o 02_trim 00_raw/{}_F.fastq.gz 00_raw/{}_R.fastq.gz &> {}_tg.log
```
```
# BarcodeIDChecker
BarcodeIDChecker.py -f1 {}_F_trim.fastq.gz -f2 {}_R_trim.fastq.gz -r 2020-09_GenbankSnakeMito.fasta -o {} -c 16 -M 55G &> BarcodeIDChecker.log
```

# Output
If MitoZ works (for original assembly or annotation after MITGARD), then you should expect the following output files. This includes a summary of the blast results (mean percent identity to different species), the raw blast results, the mitochondrial genome, MitoZ annotation results, and phylogenies for each gene. The log file (or STDOUT if log file not saved) will print each phylogeny. 
```
BarcodeIDChecker.log
BarcodeIDChecker_results/
├── blast_query.fasta
├── blast_results.summarized
├── blast_results.tab
├── {}_mitogenome.fasta
├── mitoz_most_related_species.txt
├── mitoz.result
│   ├── {}.cds
│   ├── {}.circos.dep
│   ├── {}.circos.karyotype.txt
│   ├── {}.circos.png
│   ├── {}.circos.svg
│   ├── {}.errorsummary.val
│   ├── {}.fasta
│   ├── {}.misc_feature
│   ├── {}_mitoscaf.fa.gbf
│   ├── {}_mitoscaf.fa.sqn
│   ├── {}_mitoscaf.fa.tbl
│   ├── {}_mitoscaf.fa.val
│   ├── {}.rrna
│   ├── {}.trna
│   ├── README.txt
│   ├── summary.txt
│   └── work<kmer>.files
└── [1020K]  Phylogenetics
    ├── gene.fasta
    ├── gene.fasta.aln
    ├── gene.fasta.contree
    ├── gene.fasta.iqtree
    └── gene.fasta.trim
```

If MitoZ fails, then you should expect the following output files. Noteably, a reference fasta for MITGARD will be used and additional files for the reference chosen will be found. Additionally, instead of a phylogeny for each gene you will find a single phylogeny using the full mitochondrial genome. 
```
BarcodeIDChecker.log
BarcodeIDChecker_results/
├──  alternate_references.txt 
├──  best_reference.fasta
├──  best_reference.txt
├──  blast_results.summarized
├──  blast_results.tab
├──  {}_mitogenome.fasta
└──  Phylogenetics
    ├──  fullgenome.fasta
    ├──  fullgenome.fasta.aln
    ├──  fullgenome.fasta.contree
    ├──  fullgenome.fasta.iqtree
    └──  fullgenome.fasta.trim
```

# Cite
Because this program only works as a wrapper for other programs, we recommend that you cite them as well. 
- [BarcodeIDChecker](https://github.com/reptilerhett/BarcodeIDChecker)
- [MitoZ](https://github.com/linzhi2013/MitoZ)
- [MITGARD](https://github.com/pedronachtigall/MITGARD)
  - [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
  - [Spades](https://cab.spbu.ru/software/rnaspades/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Samtools](http://www.htslib.org/)
- [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [PANDAS](https://pandas.pydata.org/)
- [dfply](https://github.com/kieferk/dfply)
- [BioPython](https://biopython.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Trimal](http://trimal.cgenomics.org/)
- [IQ-TREE](http://www.iqtree.org/)
