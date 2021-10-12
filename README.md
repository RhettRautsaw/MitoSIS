
![](logo.png)

## Rhett M. Rautsaw & Pedro G. Nachtigall 

 `MitoSIS` (Mitochondrial Species Identification System) is a wrapper for mitochondrial genome assembly and identification of sample contamination or mislabeling. Specifically, `MitoSIS` maps raw or trimmed reads to a database of reference mitochondrial sequences. It calculates the percentage of reads that map to different species using [`Kallisto`](https://pachterlab.github.io/kallisto/) to assess potential sample contamination. It then uses [`MITGARD`](https://github.com/pedronachtigall/MITGARD) to assemble and [`MitoZ`](https://github.com/linzhi2013/MitoZ) annotate the full mitochondrial genome. It then BLASTs the resulting mitogenome or barcoding genes (*e.g.*, CYTB, COX1, ND4, 16S, etc.) to check for sample mislabeling. Finally, `MitoSIS` uses a `MAFFT` and [`IQ-TREE`](http://www.iqtree.org/) to calculate phylogenetic distance to closely related species.

# Pipeline

- Map fastq reads to reference fasta using `kallisto`
	- Calculate total reads/tpm for each species in database
- Identify the best reference sequence
- Assemble the mitogenome using `MITGARD`
- Annotate mitogenome using `MitoZ`
	- Extract protein coding/barcoding genes
- Blast mitogenome or genes to reference database
	- Calculate mean percent identity for each species
- Align sequences and build phylogeny
	- Calculate mean/minimum phylogenetic distance for each species

![](MitoSIS_Flowchart.png)

# Arguments

|       flag      |   description   |
|-----------------|-----------------|
| -h, --help      | Show this help message and exit. | 
| -f1, --fastq1   | fastq read pair 1 (forward). **Default: None** |
| -f2, --fastq2   | fastq read pair 2 (reverse). **Default: None** |
| -s, --single    | single-end fastq. **Default: None** |
| -r, --reference | `genbank` *OR* `fasta+sp` database **Default: None** <br> See section below on fasta & custom databases <br> *Recommend downloading all mitochondrial data for your clade of interest <br> e.g., snakes; [Genbank Example](https://www.ncbi.nlm.nih.gov/nuccore/?term=snakes%5Bporgn%5D+AND+mitochondrion%5Bfilter%5D) <br> Send to > Complete Record > Genbank* |
| -o, --output    | Prefix for output files. **Default: 'ZZZ'** |
| -c, --cpu       | Number of threads to be used in each step. **Default: 8** |
| -M, --memory    | Max memory for Trinity (see [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity) for format). **Default: '30G'** |
| --clade         | Clade used for MitoZ. Options: 'Chordata' or 'Arthropoda'. **Default: 'Chordata'** |
| --convert       | Only perform Genbank to Fasta conversion and create a tab-delimited taxa id file     |
| --version       | Show program's version number and exit |
|<img width=200/> |<img width=500/>|

## Reference Databases
The user can download nucleotide sequences from the taxonomic group of interest from the database of [NCBI](https://www.ncbi.nlm.nih.gov/genbank/). For instance, the user can search for ["snakes\[porgn\]AND mitochondrion\[filter\]"](https://www.ncbi.nlm.nih.gov/nuccore/?term=snakes%5Bporgn%5D+AND+mitochondrion%5Bfilter%5D) and send all complete records to a Genbank formatted file. Then the GenBank fromat file is used as input in the option ```-r``` to be used as reference in MitoSIS pipeline.

The GenBank format file is converted into two files to generate a `fasta+sp` database, which is used in all steps of MitoSIS workflow. To improve the reference database by adding custom/private sequences, see the section below.

## Fasta & Custom Reference Databases
Fasta reference databases must be accompanied by a tab-delimited taxa id (`.sp`) file. We refer to this combination of files as a `fasta+sp` database. The tab-delimited taxa id (`.sp`) file must occur in the same directory as the fasta file and have the same filename with `.sp` appended (*i.e.*, ReferenceDB.fasta and ReferenceDB.fasta.sp).

If you have a Genbank database and only want to add additional or custom/private sequences, we recommend first running `--convert`. 

```
MitoSIS.py -r ReferenceDB.gb --convert
```

`--convert` will convert your Genbank file to a `fasta+sp` database without running the rest of MitoSIS. Output will be:

- `ReferenceDB.fasta`
- `ReferenceDB.fasta.sp`

With the initial `fasta+sp` database created...

Manually add your additional or custom sequences to the fasta and the identifer/taxa information to the `.sp` file. 

### `fasta+sp` Format
Each fasta sequences must have unique identifiers (similar to Genbank Accession Numbers) and those identifiers must match in the tab-delimited taxa id file. Ensure to not have descriptions in the fasta header (i.e., no spaces " " in the header, only the sequence id).

{ReferenceDB}.fasta 
```
>ID_1
ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
>ID_2
ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
```

{ReferenceDB}.fasta.sp
```
ID_1    Genus species
ID_2    Genus species
```


# Installation
**System Requirement**
- Linux

**Conda Installation**
```
# Clone this Repository
git clone https://github.com/RhettRautsaw/MitoSIS.git
cd MitoSIS
echo "export PATH=\$PATH:$PWD" >> ~/.bash_profile

# Clone MITGARD Repository 
git clone https://github.com/pedronachtigall/MITGARD.git
# Fix shebangs in MITGARD supporting scripts
sed -i '1 s/^.*$/\#\!\/usr\/bin\/env python/' MITGARD/bin/sam2msa.py
sed -i '1 s/^.*$/\#\!\/usr\/bin\/env python/' MITGARD/bin/msa2consensus.py
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
python MITGARD/install_NCBITaxa.py

# YOU'RE READY TO GO
# Check if MitoSIS.py is in your path
MitoSIS.py -h
```

# Example
Before running, we recommend testing `MitoSIS` with our [Tutorial](https://github.com/RhettRautsaw/MitoSIS/tree/master/Tutorial) dataset. 

We also recommend trimming your own data first prior to running this program. Example trimming using [Trim-Galore](https://github.com/FelixKrueger/TrimGalore) shown below. Depending on whether you are working with DNA or RNA-Seq data, you may want to change the length/quality parameters.
```
# Trimming
trim_galore --paired --phred33 --length 30 -q 20 -o 02_trim 00_raw/{}_F.fastq.gz 00_raw/{}_R.fastq.gz &> {}_tg.log
```

Below are outlines for running `MitoSIS`.
```
# MitoSIS - paired-end
MitoSIS.py -f1 {}_F_trim.fastq.gz -f2 {}_R_trim.fastq.gz -r ReferenceDB.gb -o {} -c 16 -M 55G &> MitoSIS.log

# MitoSIS - single
MitoSIS.py -s {}_merged.fastq.gz -r ReferenceDB.gb -o {} -c 16 -M 55G &> MitoSIS.log

# MitoSIS - paired-end & fasta+sp reference database
# NOTE: MitoSIS expects ReferenceDB.fasta.sp to occur in the same directory as ReferenceDB.fasta
MitoSIS.py -f1 {}_F_trim.fastq.gz -f2 {}_R_trim.fastq.gz -r ReferenceDB.fasta -o {} -c 16 -M 55G &> MitoSIS.log
```

# Output

The user can find a detailed results in the ```MitoSIS_summary_output.html``` with the potential contamination, percent identity and alignment distance across genes and all phylogenetic trees build. Moreover, during processing MitoSIS print messages at the terminal summarizing all results, which may also be used by the user to check the results. Find below an example of the printed message and a detailed information about all files generated by MitoSIS pipeline that can be used/analyzed *a posteriori* by the user.

```
MitoSIS.py -f1 CON45_R1.fq.gz -f2 CON45_R2.fq.gz -r ReferenceDB.gb -o CON45 -c 16 -M 62G

	       o O       o O       o O       o O       o O       o O
	     o | | O   o | | O   o | | O   o | | O   o | | O   o | | O
	   O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
	  O-oO | | o   O | | o   O | | o   O | | o   O | | o   O | | oO-o
	 O---o O o       O o       O o       O o       O o       O o O---o
	O-----O                                                     O-----o
	o-----O         ___  ____ _        _____ _____ _____        o-----O
	 o---O          |  \/  (_) |      /  ___|_   _/  ___|        o---O 
	  o-O           | .  . |_| |_ ___ \ `--.  | | \ `--.          o-O
	   O            | |\/| | | __/ _ \ `--. \ | |  `--. \          O
	  O-o           | |  | | | || (_) /\__/ /_| |_/\__/ /         O-O
	 O---o          \_|  |_/_|\__\___/\____/ \___/\____/         O---o
	O-----o                        v1.2                         O-----o
	o-----O                                                     o-----O
	 o---O o O       o O       o O       o O       o O       o O o---O
	  o-Oo | | O   o | | O   o | | O   o | | O   o | | O   o | | Oo-O
	   O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
	     O | | o   O | | o   O | | o   O | | o   O | | o   O | | o
	       O o       O o       O o       O o       O o       O o
	

2021-10-12 10:39:21 ::: starting MitoSIS...
	Forward Reads -> /zfs/venom/Rhett/bin/MitoSIS/Tutorial/CON45_R1.fq.gz
	Reverse Reads -> /zfs/venom/Rhett/bin/MitoSIS/Tutorial/CON45_R2.fq.gz
	Reference Database -> /zfs/venom/Rhett/bin/MitoSIS/Tutorial/ReferenceDB.gb
	Output -> /zfs/venom/Rhett/bin/MitoSIS/Tutorial/MitoSIS_results/CON45*
	Number of CPU -> 16
	Amount of memory -> 62G
	MitoZ Clade -> Chordata

2021-10-12 10:39:21 ::: Converting Genbank to Fasta :::


2021-10-12 10:39:21 ::: Converted 521 Genbank records to Fasta :::


2021-10-12 10:39:21 ::: Running kallisto index :::


2021-10-12 10:39:21 ::: Running kallisto :::


2021-10-12 10:39:21 ::: Summarizing kallisto to assess potential contamination :::

                species  read_count_sum  read_count_mean     tpm_sum  tpm_mean  read_sum_percent  read_mean_percent  tpm_sum_percent  tpm_mean_percent
    Crotalus adamanteus          3389.0           12.191  850947.699  3060.963            85.044             67.018           85.095            67.126
      Crotalus horridus           198.0            3.882   49427.400   969.165             4.969             21.343            4.943            21.253
 Agkistrodon piscivorus           398.0            2.117   99625.600   529.923             9.987             11.638            9.963            11.621

2021-10-12 10:39:21 ::: Identifying best reference sequence :::

                           sseqid  length  read_count       tpm                 species
 Selected Reference > NC_041524.1  17242      1693.75  423259.0     Crotalus adamanteus
                       MH626511.1  17242      1693.75  423259.0     Crotalus adamanteus
                       DQ523161.1  17213       199.00   49812.8  Agkistrodon piscivorus
                      NC_009768.1  17213       199.00   49812.8  Agkistrodon piscivorus
                      NC_014400.1  17260        99.00   24713.7       Crotalus horridus
                       HM641837.1  17260        99.00   24713.7       Crotalus horridus

2021-10-12 10:39:22 ::: Running MITGARD :::


2021-10-12 10:52:17 ::: Annotating MITGARD mitogenome with MitoZ :::


2021-10-12 10:59:57 ::: Moving onto BLAST :::


2021-10-12 10:59:57 ::: Running makeblastdb :::


2021-10-12 10:59:57 ::: Running BLAST :::


2021-10-12 10:59:58 ::: Summarizing Mean Percent Identity across genes :::

                    species  Mean_Percent_Identity
        Crotalus adamanteus              99.838545
 Crotalus horridus horridus              96.175000
          Crotalus horridus              89.570667
     Agkistrodon piscivorus              84.574576

2021-10-12 11:00:03 ::: Extracting BLAST matches for Phylogenetics :::


2021-10-12 11:00:03 ::: Aligning, Trimming, and Inferring Phylogeny for ND1.fasta :::

                species  minimum_distance  mean_distance
    Crotalus adamanteus          0.000000       0.000003
      Crotalus horridus          0.386271       0.424664
 Agkistrodon piscivorus          0.738328       0.738328
                                                   / NC_009768.1 Agkistrodon piscivorus
/--------------------------------------------------+                                   
|                                                  \ DQ523161.1 Agkistrodon piscivorus 
|                                                                                      
|                 /GAAZ01001454.1 Crotalus horridus                                    
+                 +                                                                    
|                 |GBKC01002148.1 Crotalus horridus                                    
|       /---------+                                                                    
|       |         |                                / HM641837.1 Crotalus horridus      
|       |         \--------------------------------+                                   
\-------+                                          \ NC_014400.1 Crotalus horridus     
        |                                                                              
        |                /JU175111.1 Crotalus adamanteus                               
        \----------------+                                                             
                         |MH626511.1 Crotalus adamanteus                               
                         +                                                             
                         |GBEX01002025.1 Crotalus adamanteus                           
                         +                                                             
                         |CON45                                                        
                         +                                                             
                         \NC_041524.1 Crotalus adamanteus                              

... MORE GENES ...

2021-10-12 11:01:10 ::: Summarizing Phylogenetic Distance across genes :::

                    species  minimum_distance  mean_distance
        Crotalus adamanteus          0.000000       0.071839
          Crotalus horridus          0.199598       0.479474
 Crotalus horridus horridus          0.241776       0.314130
     Agkistrodon piscivorus          0.461318       0.708560

2021-10-12 11:01:10 ::: Concatenating Genes and Removing Individuals with > 50% Missing :::


2021-10-12 11:01:12 ::: Running Concatenated Phylogeny :::

                                      /NC_009768.1_Agkistrodon_piscivorus      
                   /------------------+                                        
/------------------+                  \DQ523161.1_Agkistrodon_piscivorus       
|                  |                                                           
|                  \------------------------ EF669477.1_Agkistrodon_piscivorus 
+                                                                              
|                                          /HM641837.1_Crotalus_horridus       
|        /---------------------------------+                                   
|        |                                 \NC_014400.1_Crotalus_horridus      
\--------+                                                                     
         |                   /MH626511.1_Crotalus_adamanteus                   
         \-------------------+                                                 
                             |CON45_                                           
                             +                                                 
                             \NC_041524.1_Crotalus_adamanteus                  

2021-10-12 11:01:14 ::: Generating plots and HTML output :::


2021-10-12 11:01:20 ::: FINISHED :::

```


You should expect the following output files; however, sometimes the `mitoz.result` folder may be absent if the mitogenome assembly is not complete. Regardless, results include kallisto contamination results/summary (mean/sum read counts and tpm to each species), blast results/summary (mean percent identity to different species), the mitochondrial genome, `MitoZ` annotation results, and phylogenies/phylogenetic distance summary. The log file (or STDOUT if log file not saved) will print each phylogeny. 

```
MitoSIS_results/
├── alternate_references.tsv
├── best_reference.fasta
├── blast_query.fasta
├── blast_results.tsv
├── blast_summary.png
├── blast_summary.tsv
├── CON45_mitogenome.fasta
├── kallisto
│   ├── abundance.h5
│   ├── abundance.tsv
│   ├── abundance_species.tsv
│   └── run_info.json
├── kallisto_contamination.png
├── kallisto_contamination.tsv
├── MitoSIS_summary_output.html
├── mitoz.result
│   ├── CON45.cds
│   ├── CON45.circos.karyotype.txt
│   ├── CON45.circos.png
│   ├── CON45.circos.svg
│   ├── CON45.errorsummary.val
│   ├── CON45.fasta
│   ├── CON45.misc_feature
│   ├── CON45_mitoscaf.fa.gbf
│   ├── CON45_mitoscaf.fa.sqn
│   ├── CON45_mitoscaf.fa.tbl
│   ├── CON45_mitoscaf.fa.val
│   ├── CON45.rrna
│   ├── CON45.trna
│   └── summary.txt
├── phylogenetic_distance_summary.png
├── phylogenetic_distance_summary.tsv
├── Phylogenetics
│   ├── gene.fasta
│   ├── gene.fasta.aln
│   ├── gene.fasta.contree
│   ├── gene.fasta.iqtree
│   ├── gene.fasta.phylodist.tsv
│   ├── gene.fasta.png
│   ├── gene.fasta.trim
│   ├── gene.nex
│   ├── Concatenated.nex
│   ├── Concatenated.phy
│   ├── Concatenated.phy.bionj
│   ├── Concatenated.phy.ckp.gz
│   ├── Concatenated.phy.contree
│   ├── Concatenated.phy.iqtree
│   ├── Concatenated.phy.log
│   ├── Concatenated.phy.mldist
│   ├── Concatenated.phy.model.gz
│   ├── Concatenated.phy.png
│   ├── Concatenated.phy.splits.nex
│   ├── Concatenated.phy.treefile
│   └── Concatenated.phy.uniqueseq.phy
└── RearrangementCheck
    ├── align0.sam
    ├── consensus0.mfa.fasta
    ├── contigs0
    ├── newref0
    └── newref0_mitogenome.fa
```

# Cite
Because this program only works as a wrapper for other programs, we recommend that you cite them as well. 
- [MitoSIS](https://github.com/RhettRautsaw/MitoSIS)
- [PANDAS](https://pandas.pydata.org/)
- [dfply](https://github.com/kieferk/dfply)
- [Kallisto](https://pachterlab.github.io/kallisto/)
- [MITGARD](https://github.com/pedronachtigall/MITGARD)
  - [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
  - [Spades](https://cab.spbu.ru/software/rnaspades/)
- [MitoZ](https://github.com/linzhi2013/MitoZ)
- [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [BioPython](https://biopython.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Trimal](http://trimal.cgenomics.org/)
- [IQ-TREE](http://www.iqtree.org/)
- [dendropy](https://dendropy.org/)
- [NUMPY](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
