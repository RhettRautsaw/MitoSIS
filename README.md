
![](logo.png)

## Rhett M. Rautsaw & Pedro Nachtigall 

 `MitoSIS` (Mitochondrial Species Identification System) is a wrapper for mitochondrial genome assembly and identification of sample contamination or mislabeling. Specifically, `MitoSIS` maps raw or trimmed reads to a database of reference mitochondrial sequences. It calculates the percentage of reads that map to different species using [`Kallisto`](https://pachterlab.github.io/kallisto/) to assess potential sample contamination. It then uses [`MitoZ`](https://github.com/linzhi2013/MitoZ) and [`MITGARD`](https://github.com/pedronachtigall/MITGARD) to assemble and annotate the full mitochondrial genome and BLASTs the resulting mitogenome or barcoding genes (*e.g.*, CYTB, COX1, ND4, 16S, etc.) to check for sample mislabeling. Finally, `MitoSIS` uses a `MAFFT` and [`IQ-TREE`](http://www.iqtree.org/) to calculate alignment distance and infer a phylogeny.

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
| --clade         | Clade used for MitoZ. Options: 'Chordata' or 'Arthropoda'. **Default is 'Chordata'** |
| --version       | Show program's version number and exit |
|<img width=200/> |<img width=500/>|

# Pipeline

1. Map fastq reads to reference using `bwa` and `kallisto`
2. Keep reads that successfully mapped
3. Assemble mitogenome using `MitoZ`
4. If `MitoZ` fails, identify the best reference sequence
	- 4.1 Use `MITGARD` to assemble mitogenome
	- 4.2 Use `MitoZ` to annotate mitogenome (if not already done)
5. Extract protein coding/barcoding genes
6. Blast mitogenome or genes to reference database
7. Export results and sequences
8. Align sequences and build phylogeny 

![](MitoSIS_Flowchart.png)


# Installation
**System Requirement**
- Linux

**Conda Installation**
```
# Clone this Repository
git clone https://github.com/reptilerhett/MitoSIS.git
cd MitoSIS
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
conda env create -f mitosis_env.yml
conda activate mitosis_env

# Install dfply
pip install dfply

# Install Taxonomy Database for MitoZ
python
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
quit()

# YOU'RE READY TO GO
# Check if MitoSIS.py is in your path
MitoSIS.py -h
```

# Example
We recommend trimming your data first prior to running this program. Example trimming using [Trim-Galore](https://github.com/FelixKrueger/TrimGalore) shown below. Depending on whether you are working with DNA or RNA-Seq data, you may want to change the length/quality parameters.
```
# Trimming
trim_galore --paired --phred33 --length 30 -q 20 -o 02_trim 00_raw/{}_F.fastq.gz 00_raw/{}_R.fastq.gz &> {}_tg.log
```
```
# MitoSIS
MitoSIS.py -f1 {}_F_trim.fastq.gz -f2 {}_R_trim.fastq.gz -r 2020-09_GenbankSnakeMito.gb -o {} -c 16 -M 55G &> MitoSIS.log
```

# Output
```
       o O       o O       o O       o O       o O       o O
     o | | O   o | | O   o | | O   o | | O   o | | O   o | | O
   O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
  O-oO | | o   O | | o   O | | o   O | | o   O | | o   O | | oO-o
 O---o O o       O o       O o       O o       O o       O o O---o
O-----O                                                     O-----o
o-----O                                                     o-----O
 o---O   ,--.   ,--.,--.  ,--.          ,---.  ,--. ,---.    o---O 
  o-O    |   `.'   |`--',-'  '-. ,---. '   .-' |  |'   .-'    o-O
   O     |  |'.'|  |,--.'-.  .-'| .-. |`.  `-. |  |`.  `-.     O
  O-o    |  |   |  ||  |  |  |  ' '-' '.-'    ||  |.-'    |   O-O
 O---o   `--'   `--'`--'  `--'   `---' `-----' `--'`-----'   O---o
O-----o                                                     O-----o
o-----O                                                     o-----O
 o---O o O       o O       o O       o O       o O       o O o---O
  o-Oo | | O   o | | O   o | | O   o | | O   o | | O   o | | Oo-O
   O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
     O | | o   O | | o   O | | o   O | | o   O | | o   O | | o
       O o       O o       O o       O o       O o       O o


2020-11-21 17:36:46 >>>> starting MitoSIS...
	Forward Reads -> /zfs/venom/Rhett/2020_BarcodeTest/CON01/CON01_R1.fq
	Reverse Reads -> /zfs/venom/Rhett/2020_BarcodeTest/CON01/CON01_R2.fq
	Reference Database -> /zfs/venom/Rhett/2020_BarcodeTest/2020-10_GenbankSnakeMito/2020-10_SnakeMito.fasta
	Output -> /zfs/venom/Rhett/2020_BarcodeTest/CON01/MitoSIS_results/CON01*
	Number of CPU -> 16
	Amount of memory -> 55G


2020-11-21 17:36:46 ::: Genbank to Fasta conversion previously completed :::


2020-11-21 17:36:48 ::: kallisto index previously completed :::


2020-11-21 17:36:48 ::: Running kallisto :::


2020-11-21 17:36:53 ::: Summarizing kallisto to assess potential contamination :::

             species  read_count  read_percent
 Crotalus adamanteus  238.000046     98.755187
   Crotalus horridus    3.000000      1.244813

2020-11-21 17:37:00 ::: bwa index previously completed :::


2020-11-21 17:37:00 ::: Running bwa mem :::


2020-11-21 17:37:01 ::: Sorting/converting sam files :::


2020-11-21 17:37:01 ::: Converting bam to fastq :::


2020-11-21 17:37:01 ::: Running MitoZ assembly :::


2020-11-21 17:37:02 ::: MitoZ failed, trying MITGARD :::


2020-11-21 17:37:02 ::: Identifying best reference sequence :::

      sseqid  length  read_count              species
  MH626511.1   17242       119.0  Crotalus adamanteus
 NC_041524.1   17242       119.0  Crotalus adamanteus
 NC_014400.1   17260         1.5    Crotalus horridus
  HM641837.1   17260         1.5    Crotalus horridus

2020-11-21 17:37:03 ::: Running MITGARD :::


2020-11-21 17:38:21 ::: Annotating MITGARD mitogenome with MitoZ :::


2020-11-21 17:46:02 ::: Moving onto BLAST :::


2020-11-21 17:46:02 ::: makeblastdb previously completed :::


2020-11-21 17:46:02 ::: Running blast :::


2020-11-21 17:46:02 ::: Summarizing Mean Percent Identity across genes :::

                     species  Mean_Percent_Identity
         Crotalus adamanteus              99.808421
 Crotalus oreganus caliginis              90.370000
         Crotalus mitchellii              90.080000
   Crotalus oreganus abyssus              90.040000
   Crotalus oreganus helleri              89.746667
           Crotalus cerberus              89.698511
    Crotalus viridis nuntius              89.690000

2020-11-21 17:46:06 ::: Extracting BLAST matches for Phylogenetics :::


2020-11-21 17:46:06 ::: Aligning, Trimming, and Inferring Phylogeny for ND1.fasta :::

                             , GBEX01002025.1_TSA__Crotalus_adamante...
                             |
                             | CON01_CON01_ND1_len_981__2537_3517____
                             |
                             | NC_041524.1_Crotalus_adamanteus_mitoc...
                             |
                  ___________| MH626511.1_Crotalus_adamanteus_mitoch...
                 |           |
                 |           | JU175111.1_TSA__Crotalus_adamanteus_C...
             ____|
            |    |               , NC_014400.1_Crotalus_horridus_mitocho...
            |    |            ___|
            |    |           |   | HM641837.1_Crotalus_horridus_mitochon...
  __________|    |___________|
            |                , GBKC01002148.1_TSA__Crotalus_horridus...
            |                |
            |                | GAAZ01001454.1_TSA__Crotalus_horridus...
            |
            |__________________ MK313588.1_Crotalus_cerastes_voucher_...


2020-11-21 17:50:01 ::: Summarizing Mean Alignment Distance across genes :::

                     species  Mean_Alignment_Distance
         Crotalus adamanteus                 0.001493
           Crotalus cerberus                 0.102720
           Crotalus horridus                 0.104204
 Crotalus oreganus caliginis                 0.105366
   Crotalus oreganus helleri                 0.105685
         Crotalus mitchellii                 0.108293
  Crotalus oreganus oreganus                 0.111545

2020-11-21 17:50:01 ::: Concatenating Genes and Removing Individuals with > 50% Missing :::


2020-11-21 17:50:04 ::: Running Concentated Phylogeny :::

                                             , CON01_
  ___________________________________________|
 |                                           , NC_041524.1_Crotalus_adamanteus
 |                                           |
_|                                           | MH626511.1_Crotalus_adamanteus
 |
 |                                           , NC_014400.1_Crotalus_horridus
 |___________________________________________|
                                             | HM641837.1_Crotalus_horridus


2020-11-21 17:50:05 ::: FINISHED :::
```

If `MitoZ` is successful (for original assembly or annotation after `MITGARD`), then you should expect the following output files. This includes a summary of the blast results (mean percent identity to different species), the raw blast results, the mitochondrial genome, `MitoZ` annotation results, and phylogenies for each gene. The log file (or STDOUT if log file not saved) will print each phylogeny. 
```
MitoSIS.log
MitoSIS_results/
├── alignment_summary.tsv
├── blast_query.fasta
├── blast_summary.tsv
├── blast_results.tsv
├── kallisto_contamination.tsv
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
└── Phylogenetics
    ├── gene.fasta
    ├── gene.fasta.aln
    ├── gene.fasta.contree
    ├── gene.fasta.iqtree
    └── gene.fasta.trim
```

If `MitoZ` fails, then you should expect the following output files. Noteably, a reference mitogenome will be chosen and export. Additionally, instead of a phylogeny for each gene you will find a single phylogeny using the full mitochondrial genome. 
```
MitoSIS.log
MitoSIS_results/
├──  alignment_summary.tsv
├──  alternate_references.tsv 
├──  best_reference.fasta
├──  best_reference.txt
├──  blast_query.fasta
├──  blast_summary.tsv
├──  blast_results.tsv
├──  kallisto_contamination.tsv
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
- [MitoSIS](https://github.com/reptilerhett/MitoSIS)
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
- [Kallisto](https://pachterlab.github.io/kallisto/)
