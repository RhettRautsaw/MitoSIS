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
