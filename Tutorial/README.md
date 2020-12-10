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

The user can check the expected output printed at the terminal in the following files: [`ExpectedOut_paired.md`](https://github.com/reptilerhett/MitoSIS/blob/master/Tutorial/ExpectedOut_paired.md) [`ExpectedOut_single.md`](https://github.com/reptilerhett/MitoSIS/blob/master/Tutorial/ExpectedOut_single.md).

The printed messages were designed to show a summary of all results obtained with MitoSIS.

Moreover, MitoSIS compiles all results obtained in an user-friendly HTML file ```MitoSIS_summary_output.html``` containing charts and the phylogenetic trees. Also, the user can use any file generated during MitoSIS processing in any downstream analysis by checking the ```MitoSIS_output``` directory.
