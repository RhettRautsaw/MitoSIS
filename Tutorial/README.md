# Tutorial

A quick tutorial to use `MitoSIS`. 

Data provided in this repository (`CON45_*.fq.gz`) is 4,000 paired-end reads simulated from `NC_041524` (*Crotalus adamanteus*) with 5% contamination from `NC_014400` (*Crotalus horridus*) and 10% contamination from `NC_009768` (*Agkistrodon piscivorus*).

The `ReferenceDB.gb` is a reference database with all mitochondrial sequences from *Crotalus adamanteus*, *Crotalus horridus*, and *Agkistrodon piscivous*. You will first need to gunzip this file.

```
gunzip ReferenceDB.gb.gz
```

# Converting `Genbank` to `fasta+sp`
This step is done automatically in MitoSIS when you provide a Genbank file, but you can also do the conversion independent of the rest of `MitoSIS`.

```
MitoSIS.py -r ReferenceDB.gb --convert
```

# Running Paired-End

```
MitoSIS.py -f1 CON45_R1.fq.gz -f2 CON45_R2.fq.gz -r ReferenceDB.gb -o CON45 -c 8 -M 16G --clade Chordata
```

# Running Single-End
If the user wants to use/test single-end mode, reads can concatenated together. Reads could also be merged using [`PEAR`](https://cme.h-its.org/exelixis/web/software/pear/). For this tutorial, we choose to break the pairs and concatenate files together into one file. 
```
cat CON45_R1.fq.gz CON45_R2.fq.gz > CON45_merge.fq.gz

MitoSIS.py -s CON45_merge.fq.gz -r ReferenceDB.gb -o CON45 -c 8 -M 16G --clade Chordata
```

# Expected Output

The user can check the expected output by decompressing the `.tar.gz` files. 
```
tar -xvzf *.tar.gz
```

MitoSIS compiles all results obtained in an user-friendly HTML file `MitoSIS_summary_output.html` containing charts and the phylogenetic trees. Also, the user can use any file generated during MitoSIS processing in any downstream analysis by checking the `MitoSIS_output` directory.

Morover, the printed messages were designed to show a summary of all results obtained with MitoSIS. Here's the output for the paired-end tutorial...

```

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
	O-----o                                                     O-----o
	o-----O                                                     o-----O
	 o---O o O       o O       o O       o O       o O       o O o---O
	  o-Oo | | O   o | | O   o | | O   o | | O   o | | O   o | | Oo-O
	   O | | | | O | | | | O | | | | O | | | | O | | | | O | | | | O
	     O | | o   O | | o   O | | o   O | | o   O | | o   O | | o
	       O o       O o       O o       O o       O o       O o
	

2020-12-23 12:09:01 ::: starting MitoSIS...
	Forward Reads -> /scratch1/rrautsa/00_PedroMitoSISTest/MitoSIS/Tutorial/CON45_R1.fq.gz
	Reverse Reads -> /scratch1/rrautsa/00_PedroMitoSISTest/MitoSIS/Tutorial/CON45_R2.fq.gz
	Reference Database -> /scratch1/rrautsa/00_PedroMitoSISTest/MitoSIS/Tutorial/ReferenceDB.gb
	Output -> /scratch1/rrautsa/00_PedroMitoSISTest/MitoSIS/Tutorial/MitoSIS_results/CON45*
	Number of CPU -> 16
	Amount of memory -> 62G
	MitoZ Clade -> Chordata


2020-12-23 12:09:01 ::: Genbank to Fasta conversion previously completed :::


2020-12-23 12:09:01 ::: kallisto index previously completed :::


2020-12-23 12:09:01 ::: Running kallisto :::


2020-12-23 12:09:01 ::: Summarizing kallisto to assess potential contamination :::

                species  read_count         tpm  read_percent  tpm_percent
    Crotalus adamanteus      3389.0  850947.699        85.044       85.095
 Agkistrodon piscivorus       398.0   99625.600         9.987        9.963
      Crotalus horridus       198.0   49427.400         4.969        4.943

2020-12-23 12:09:01 ::: bwa index previously completed :::


2020-12-23 12:09:01 ::: Running bwa mem :::


2020-12-23 12:09:01 ::: Sorting/converting sam files :::


2020-12-23 12:09:04 ::: Converting bam to fastq :::


2020-12-23 12:09:14 ::: Running MitoZ assembly :::


2020-12-23 12:24:10 ::: MitoZ ran successfully :::


2020-12-23 12:24:10 ::: Moving onto BLAST :::


2020-12-23 12:24:10 ::: makeblastdb previously completed :::


2020-12-23 12:24:10 ::: Running BLAST :::


2020-12-23 12:24:10 ::: Summarizing Mean Percent Identity across genes :::

                    species  Mean_Percent_Identity
        Crotalus adamanteus              99.838545
 Crotalus horridus horridus              96.175000
          Crotalus horridus              89.570667
     Agkistrodon piscivorus              84.574576

2020-12-23 12:24:10 ::: Extracting BLAST matches for Phylogenetics :::


2020-12-23 12:24:10 ::: Aligning, Trimming, and Inferring Phylogeny for ND1.fasta :::

                                   , NC_009768.1_Agkistrodon_piscivorus_mi...
  _________________________________|
 |                                 | DQ523161.1_Agkistrodon_piscivorus_mit...
 |
 |                            , NC_041524.1_Crotalus_adamanteus_mitoc...
 |                            |
 |                            | CON45_CON45_ND1_len_981__2536_3516____
_|                            |
 |                            | GBEX01002025.1_TSA__Crotalus_adamante...
 |                            |
 |            ________________| MH626511.1_Crotalus_adamanteus_mitoch...
 |           |                |
 |           |                | JU175111.1_TSA__Crotalus_adamanteus_C...
 |___________|
             |                     , NC_014400.1_Crotalus_horridus_mitocho...
             |                _____|
             |               |     | HM641837.1_Crotalus_horridus_mitochon...
             |_______________|
                             , GBKC01002148.1_TSA__Crotalus_horridus...
                             |
                             | GAAZ01001454.1_TSA__Crotalus_horridus...

< ... removed output from other genes ... >

2020-12-23 12:25:25 ::: Summarizing Mean Alignment Distance across genes :::

                    species  Mean_Alignment_Distance
        Crotalus adamanteus                 0.023411
          Crotalus horridus                 0.162659
     Agkistrodon piscivorus                 0.200809
 Crotalus horridus horridus                 0.399161

2020-12-23 12:25:25 ::: Concatenating Genes and Removing Individuals with > 50% Missing :::


2020-12-23 12:25:27 ::: Running Concatenated Phylogeny :::

                                      _____ EF669477.1_Agkistrodon_piscivorus
  ___________________________________|
 |                                   |   , DQ523161.1_Agkistrodon_piscivorus
 |                                   |___|
 |                                       | NC_009768.1_Agkistrodon_piscivorus
_|
 |                                     , NC_041524.1_Crotalus_adamanteus
 |                                     |
 |               ______________________| CON45_
 |              |                      |
 |______________|                      | MH626511.1_Crotalus_adamanteus
                |
                |                         , NC_014400.1_Crotalus_horridus
                |_________________________|
                                          | HM641837.1_Crotalus_horridus


2020-12-23 12:25:29 ::: Generating plots and HTML output :::


2020-12-23 12:25:36 ::: FINISHED :::

```
