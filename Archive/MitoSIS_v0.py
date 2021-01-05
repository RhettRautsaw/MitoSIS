#!/usr/bin/env python

import sys, os, shutil
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio import Phylo
from Bio import AlignIO
from Bio.Nexus import Nexus
from Bio.Alphabet import IUPAC
from Bio.Phylo.TreeConstruction import DistanceCalculator
calculator=DistanceCalculator("identity")
import pandas as pd
from dfply import *
import glob

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
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

MitoSIS is a wrapper for mitochondrial genome assembly and identification of sample contamination or mislabeling. Specifically MitoSIS maps raw or trimmed reads to a database of reference mitochondrial sequences. It calculates the percentage of reads that map to different species using Kallisto to assess potential sample contamination. It then uses MitoZ and MITGARD to assemble and annotate the full mitochondrial genome and BLASTs the resulting mitogenome or barcoding genes (e.g., CYTB, COX1, ND4, 16S, etc.) to check for sample mislabeling. Finally, MitoSIS uses a MAFFT and IQTREE to calculate alignment distance and infer a phylogeny.

:: PIPELINE ::
1. Map fastq reads to reference fasta using bwa and kallisto
2. Keep reads that successfully mapped
3. Assemble mitogenome using MitoZ
4. If MitoZ fails, identify the best reference sequence
	4.1 Use MITGARD to assemble mitogenome
	4.2 Use MitoZ to annotate mitogenome (if not already done)
5. Extract protein coding/barcoding genes
6. Blast mitogenome or genes to reference database
7. Export results and sequences
8. Align sequences and build phylogeny 

:: EXAMPLE ::
MitoSIS.py -f1 sample_F.fastq.gz -f2 sample_R.fastq.gz -r 2020-09_GenbankSnakeMito.fasta -o sample -c 16 -M 55G

:: CITE :: 
https://github.com/RhettRautsaw/MitoSIS\n\n""")

###############################################

parser.add_argument("-f1","--fastq1",
					type=argparse.FileType('r+'),
					help="REQUIRED: Full path to fastq read pair 1 (forward)")
parser.add_argument("-f2","--fastq2",
					type=argparse.FileType('r+'),
					help="REQUIRED: Full path to fastq read pair 2 (reverse)")
parser.add_argument("-r","--reference",
					type=argparse.FileType('r+'),
					help="REQUIRED: Full path to fasta reference database (I recommend downloading all mitochondrial data for your clade of interest from GenBank, e.g. snakes)")
parser.add_argument("-o","--output",
					type=str,
					default='ZZZ',
					help="OPTIONAL: Prefix for output files. Default is 'ZZZ'")
parser.add_argument("-c","--cpu",
					type=int,
					default=8,
					help="OPTIONAL: Number of threads to be used in each step. Default is 8")
parser.add_argument("-M","--memory",
					type=str,
					default='30G',
					help="OPTIONAL: Max memory for Trinity assembler, use the same format as Trinity. Default is '30G'")
parser.add_argument("--clade",
					type=str,
					default='Chordata',
					help="Clade used for MitoZ. Options: 'Chordata' or 'Arthropoda'. Default is 'Chordata'")
parser.add_argument("--version", action='version', version='MitoSIS v1')
args=parser.parse_args()

############################################### SETUP

#fastq1_name="/zfs/venom/Rhett/Data/SeqCap/I0771_Agkistrodon_conanti/02_trim/I0771_Agkistrodon_conanti_F_trim.fastq.gz"
#fastq2_name="/zfs/venom/Rhett/Data/SeqCap/I0771_Agkistrodon_conanti/02_trim/I0771_Agkistrodon_conanti_R_trim.fastq.gz"
#output="I0771_Agkistrodon_conanti"

#fastq1_name="/zfs/venom/Rhett/Data/SeqCap/I0796_Daboia_russelii/00_raw/I0796_Daboia_russelii_F.fastq.gz"
#fastq2_name="/zfs/venom/Rhett/Data/SeqCap/I0796_Daboia_russelii/00_raw/I0796_Daboia_russelii_R.fastq.gz"
#fastq1_name="/zfs/venom/Rhett/Data/SeqCap/I0796_Daboia_russelii/02_trim/I0796_Daboia_russelii_F_trim.fastq.gz"
#fastq2_name="/zfs/venom/Rhett/Data/SeqCap/I0796_Daboia_russelii/02_trim/I0796_Daboia_russelii_R_trim.fastq.gz"
#output="I0796_Daboia_russelii"
#num_threads=16
#memory="55G"
#
#fastq1_name="/zfs/venom/Rhett/2020_BarcodeTest/test2/CON01_R1.fq"
#fastq2_name="/zfs/venom/Rhett/2020_BarcodeTest/test2/CON01_R2.fq"
#output="CON01"
#num_threads=16
#memory="55G"
#reference_name_gb="/zfs/venom/Rhett/2020_BarcodeTest/2020-10_GenbankSnakeMito/2020-10_SnakeMito.gb"
#reference_name = reference_name_gb.split(".gb")[0]+".fasta"
#references = list(SeqIO.parse(reference_name,"fasta"))
#
#reference_name_gb="/zfs/venom/Rhett/Data/2020-10_GenbankSnakeMito/2020-10_SnakeMito.gb"
#reference_name_gb="/zfs/venom/Rhett/Data/SeqCap/I0771_Agkistrodon_conanti/bc2/2020-10_SnakeMito.gb"
#reference_name = reference_name_gb.split(".gb")[0]+".fasta"
#references = list(SeqIO.parse(reference_name,"fasta"))

fastq1_name = os.path.abspath(args.fastq1.name)
fastq2_name = os.path.abspath(args.fastq2.name)
reference_name_gb = os.path.abspath(args.reference.name)
reference_name = reference_name_gb.split(".gb")[0]+".fasta"
output = args.output
num_threads = args.cpu
memory=args.memory
clade=args.clade

print("""
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
""")
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting MitoSIS...")
CWD = os.getcwd()
print("\tForward Reads -> "+fastq1_name)
print("\tReverse Reads -> "+fastq2_name)
print("\tReference Database -> "+reference_name_gb)
print("\tOutput -> " + CWD + "/MitoSIS_results/"+output+"*")
print("\tNumber of CPU -> "+str(num_threads))
print("\tAmount of memory -> "+memory+"\n")
print("\tMitoZ Clade -> "+clade+"\n")

os.mkdir("MitoSIS_results")
os.chdir('MitoSIS_results')

############################################### PARSE GENBANK

if os.path.isfile(reference_name + ".sp"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Genbank to Fasta conversion previously completed :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting Genbank to Fasta :::\n")
	sp.call('rm ' + reference_name, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	species=[]
	species2=[]
	count=1
	gb = open(reference_name_gb, 'r')
	for gb_record in SeqIO.parse(gb,'genbank'):
		acc = gb_record.id
		org = gb_record.annotations['organism']
		acc_org = [acc,org]
		species.append(acc_org)
		for feature in gb_record.features:
			if 'db_xref' in feature.qualifiers:
				taxid=', '.join(feature.qualifiers['db_xref'])
				if 'taxon' in taxid:
					taxid=re.sub(".*taxon\:","", taxid)
					acc_taxid=[acc,taxid]
					species2.append(acc_taxid)
		with open(reference_name, "a") as output_handle:
			tmp=SeqIO.write(gb_record, output_handle, "fasta")
			count = 1 + count
	
	species_df = pd.DataFrame(species)
	species_df.to_csv(reference_name+'.sp', index=False, header=False, sep="\t")
	
	species_df2 = pd.DataFrame(species2)
	species_df2 = species_df2.drop_duplicates()
	species_df2.to_csv(reference_name+'.sp2', index=False, header=False, sep="\t")
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converted %i Genbank records to Fasta :::\n" % count)

references = list(SeqIO.parse(reference_name,"fasta"))
species = pd.read_csv(reference_name+'.sp',sep="\t", names=['sseqid','species']) 

############################################### RUN KALLISTO

if os.path.isfile(reference_name + ".kallisto"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: kallisto index previously completed :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running kallisto index :::\n")
	sp.call('kallisto index -i ' + reference_name + ".kallisto " + reference_name, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running kallisto :::\n")
sp.call("kallisto quant -i " + reference_name + ".kallisto -o kallisto -t " + str(num_threads) + " " + fastq1_name + " " + fastq2_name, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Summarizing kallisto to assess potential contamination :::\n")
kallisto = pd.read_csv("kallisto/abundance.tsv", sep="\t", names=['sseqid','length','eff_length','read_count','tpm'], skiprows=1)
kallisto = pd.merge(kallisto, species, on ='sseqid', how ='left')
kallisto_sum=kallisto >> group_by(X.species) >> summarize(read_count = X.read_count.sum(), tpm = X.tpm.sum()) >> mask(X.read_count > 0) >> ungroup() >> mutate(read_percent=(X.read_count/X.read_count.sum())*100, tpm_percent=(X.tpm/X.tpm.sum())*100)
kallisto_sum=kallisto_sum.sort_values("tpm_percent", ascending=False)
kallisto_sum=kallisto_sum.round(3)
kallisto_sum2=kallisto_sum >> mask(X.tpm_percent > 1)
print(kallisto_sum2.to_string(index=False))
kallisto_sum.to_csv("kallisto_contamination.tsv",sep='\t',index=False)

############################################### RUN BWA

if os.path.isfile(reference_name + ".amb"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: bwa index previously completed :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running bwa index :::\n")
	sp.call('bwa index ' + reference_name, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running bwa mem :::\n")
sp.call("bwa mem -t " + str(num_threads) + " " + reference_name + " " + fastq1_name + " " + fastq2_name + " > tmp.sam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### EXTRACT MAPPED READS WITH SAMTOOLS

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Sorting/converting sam files :::\n")
sp.call("samtools sort -@ "+str(num_threads)+" tmp.sam > tmp.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

#print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Pulling out reads that successfully mapped  :::\n")
#sp.call("samtools view -b -F 4 tmp.bam > tmp2.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting bam to fastq :::\n")
sp.call("samtools fastq -@ "+str(num_threads)+" -F 4 -1 tmp_F.fq -2 tmp_R.fq -n tmp.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### RUN MITOZ

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running MitoZ assembly :::\n")
sp.call("MitoZ.py all2 --genetic_code auto --clade "+clade+" --outprefix " + output + " --thread_number " + str(num_threads) + " --fastq1 tmp_F.fq --fastq2 tmp_R.fq --run_mode 2", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### IF MITOZ FAILED, RUN MITGARD

if os.path.isdir(output+".result"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: MitoZ ran successfully :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: MitoZ failed, trying MITGARD :::\n")
	sp.call("rm -rf " + output + ".tmp/", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	############################################### IDENTIFY BEST REFERENCE SEQUENCE BY NUMBER OF MAPPED READS
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Identifying best reference sequence :::\n")
	
	alt_ref=kallisto.sort_values("read_count", ascending=False) >> mask(X.length > 10000, X.read_count > 0) >> drop(X.eff_length, X.tpm)
	
	if alt_ref.shape[0]==0:
		sp.call("samtools index tmp2.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call("samtools idxstats tmp2.bam | sort -k3nr > alternate_references.tsv", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		alt_ref = pd.read_csv("alternate_references.tsv",sep='\t',names=["sseqid", "length", "read_count", "unmapped"])
		alt_ref = pd.merge(alt_ref, species, on ='sseqid', how ='left')
		alt_ref=alt_ref >> mask(X.length > 10000, X.read_count>0) >> drop(X.unmapped)
		best_ref=alt_ref['sseqid'].iloc[0]
		print(alt_ref.to_string(index=False))
		alt_ref.to_csv("alternate_references.tsv",sep='\t',index=False)
	else:
		best_ref=alt_ref['sseqid'].iloc[0]
		alt_ref=alt_ref.round(3)
		alt_ref2=alt_ref.head(7)
		alt_ref2=alt_ref2.reset_index(drop=True)
		alt_ref2.loc[0,'sseqid']="Selected Reference > "+alt_ref2.loc[0,'sseqid']
		print(alt_ref2.to_string(index=False))
		alt_ref.to_csv("alternate_references.tsv",sep='\t',index=False)
		
	ref_seq = []
	for seq in references:
		if seq.id == best_ref:
			ref_seq.append(seq)
	
	handle=open('best_reference.fasta', "w")
	writer = FastaIO.FastaWriter(handle)
	writer.write_file(ref_seq)
	handle.close()
	
	sp.call("sed -i 's/ .*//g' best_reference.fasta", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	############################################### RUN MITGARD
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running MITGARD :::\n")
	sp.call("MITGARD.py -s tmp -1 " + fastq1_name + " -2 " + fastq2_name + " -R best_reference.fasta -r True -c " + str(num_threads) + " -M " + memory, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	sp.call("sed -i 's/>.*/>scaffold1/g' tmp_mitogenome.fa", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	############################################### ANNOTATE MITGARD ASSEMBLY WITH MITOZ
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Annotating MITGARD mitogenome with MitoZ :::\n")
	sp.call("MitoZ.py annotate --genetic_code auto --clade "+clade+" --outprefix " + output + " --thread_number " + str(num_threads) + " --fastafile tmp_mitogenome.fa", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### CHECK IF MITOZ RAN AND MOVE FILES

if os.path.isdir(output + ".result"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Moving onto BLAST :::\n")
	sp.call('mv '+output+'.result mitoz.result', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('cat mitoz.result/' + output + '.cds mitoz.result/' + output + '.rrna > blast_query.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("perl -pi -e 's/>.*?;/>"+output+";/g' blast_query.fasta", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('cp mitoz.result/*most_related_species.txt mitoz_most_related_species.txt', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('cp mitoz.result/'+output+'.fasta '+output+'_mitogenome.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("perl -pi -e 's/>.*?;/>"+output+";/g' "+output+"_mitogenome.fasta", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: MitoZ annotation failed :::\n")
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Genome assembly too fragmented :::\n")
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: We will blast what MITGARD was able to assemble :::\n")
	sp.call("sed -i 's/>scaffold1/>"+output+"/g' tmp_mitogenome.fa", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('mv tmp_mitogenome.fa '+output+'_mitogenome.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### CREATE BLAST DATABASE FROM REFERENCE SEQUENCES

if os.path.isfile(reference_name + ".nin"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: makeblastdb previously completed :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running makeblastdb :::\n")
	sp.call('makeblastdb -in ' + reference_name + ' -dbtype nucl', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### BLAST ANNOTATED REGIONS OR FULL MITOGENOME (IF ANNOTATION FAILS) TO REFERENCE DATABASE

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running blast :::\n")
if os.path.isdir("mitoz.result"):
	sp.call('blastn -query blast_query.fasta -db ' + reference_name + ' -outfmt "6 qseqid sseqid stitle pident evalue bitscore sseq" -num_threads ' + str(num_threads) + ' -max_target_seqs 50 -qcov_hsp_perc 90 -evalue 0.0001 -out blast_results.tsv', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	query = list(SeqIO.parse("blast_query.fasta","fasta"))
else:
	sp.call('blastn -query '+output+'_mitogenome.fasta -db ' + reference_name + ' -outfmt "6 qseqid sseqid stitle pident evalue bitscore sseq" -num_threads ' + str(num_threads) + ' -max_target_seqs 50 -evalue 0.0001 -out blast_results.tsv', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	query = list(SeqIO.parse(output+'_mitogenome.fasta',"fasta"))

############################################### EXTRACT BLAST MATCH SEQUENCES

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Summarizing Mean Percent Identity across genes :::\n")
header_list = ["qseqid","sseqid","stitle","pident","evalue","bitscore","sseq"]
results = pd.read_csv("blast_results.tsv",sep='\t',names=header_list)
results = pd.merge(results, species, on ='sseqid', how ='left')

results2=results >> group_by(X.species) >> summarize(Mean_Percent_Identity = X.pident.mean())
results2=results2.sort_values("Mean_Percent_Identity", ascending=False)
results3=results2.head(7)
print(results3.to_string(index=False))
results2.to_csv("blast_summary.tsv",sep='\t',index=False)

sp.call("sort -t$'\t' -k4nr blast_results.tsv > tmp.tab", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call("mv tmp.tab blast_results.tsv", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call('rm -rf kallisto tmp* *.tmp mapped bowtie_index align.sam consensus.mfa.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

os.mkdir("Phylogenetics")
os.chdir("Phylogenetics")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Extracting BLAST matches for Phylogenetics :::\n")
if os.path.isdir("../mitoz.result"):
	files=[]
	for i in results.qseqid.unique():
		gene=i.split(";")[1]
		files.append(gene+'.fasta')
		
		############################################### CONVERTING PANDAS DATAFRAME (BLAST OUTFMT 6) TO FASTA
		
		ref_seqs=results >> mask(X.qseqid==i) >> select(X.stitle, X.sseq) >> mutate(stitle=">"+X.stitle)
		fna=[]
		for index in ref_seqs.index:
			fna.append(ref_seqs['stitle'][index].replace(" "," "))
			fna.append(ref_seqs['sseq'][index].replace("-",""))
		
		with open(gene+'.fasta', 'w') as f:
			for item in fna:
				bytes=f.write('%s\n' % item)
		
		############################################### COMBINING NEW FASTA WITH QUERY SEQUENCE
		
		ref_seqs = list(SeqIO.parse(gene+'.fasta',"fasta"))
		
		gene_db = []
		for seq in query:
			if seq.id == i:
				seq.id=seq.id.split(";")[0]
				gene_db.append(seq)
		for seq in ref_seqs:
			gene_db.append(seq)
		
		############################################### WRITING COMBINED FILE TO GENE.FASTA
		
		handle=open(gene + '.fasta', "w")
		writer = FastaIO.FastaWriter(handle)
		writer.write_file(gene_db)
		handle.close()
else:
	files=[]
	files.append('fullgenome.fasta')
	gene_db=[]
	for seq in query:
		gene_db.append(seq)
	for seq in references:
		if seq.id in results.sseqid.unique():
			if len(seq.seq) > 10000:
				gene_db.append(seq)
	handle=open('fullgenome.fasta', "w")
	writer = FastaIO.FastaWriter(handle)
	writer.write_file(gene_db)
	handle.close()

############################################### ALIGNING AND MAKING A PHYLOGENY FOR EACH GENE

dist_summary=pd.DataFrame()

for i in files:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Aligning, Trimming, and Inferring Phylogeny for " + i + " :::\n")
	gene=i.split(".")[0]
	sp.call("mafft "+i+" > "+i+".aln", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("trimal -in "+i+".aln -out "+i+".trim -automated1 -keepheader", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("trimal -in "+i+".aln -out "+i+".trim2 -automated1", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("sed -i 's/ /_/g' "+i+".trim", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	AlignIO.convert(i+".trim2", "fasta", gene+".nex",'nexus',alphabet=IUPAC.ambiguous_dna)
	
	############################################### CALCULATING ALIGNMENT DISTANCE
	
	aln = AlignIO.read(open(i+'.trim2'),'fasta')
	dm = calculator.get_distance(aln)
	dist = []
	for rec in dm.matrix:
		dist.append(rec[0])
	
	df = pd.DataFrame({'sseqid' : dm.names[1:], 'dist' : dist[1:]})
	df2 = pd.merge(df, species, on ='sseqid', how ='left')
	dist_summary=dist_summary.append(df2)
	
	############################################### RUNNING IQTREE FOR EACH GENE & PRINT RESULTS
	
	sp.call("iqtree -s "+i+".trim -bb 1000 -seed 12345", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	if os.path.isfile(i+".trim.contree"):
		tree = Phylo.read(i+".trim.contree", 'newick')
		tree.root_at_midpoint()
		Phylo.write(tree,i+".contree","newick")
		Phylo.draw_ascii(tree)
		sp.call('mv '+i+'.trim.iqtree '+i+'.iqtree',shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

sp.call('rm *.trim.*', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### EXPORTING MEAN ALIGNMENT DISTANCES

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Summarizing Mean Alignment Distance across genes :::\n")
dist_summary2 = dist_summary >> group_by(X.species) >> summarize(Mean_Alignment_Distance=X.dist.mean())
dist_summary2=dist_summary2.sort_values("Mean_Alignment_Distance")
dist_summary3=dist_summary2.head(7)
print(dist_summary3.to_string(index=False))
dist_summary2.to_csv("../alignment_summary.tsv",sep='\t',index=False)

############################################### CONCATENATING GENES

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Concatenating Genes and Removing Individuals with > 50% Missing :::\n")
file_list=sorted(glob.glob('*.nex'))
nexi =  [(fname.replace("-","_"), Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('Concatenated.nex', 'w'))
combined 

############################################### REMOVING SEQUENCES WITH > 50% MISSING, GAPS, OR 'N'

sequences = list(SeqIO.parse('Concatenated.nex',"nexus"))
sequences2 = []
for seq in sequences:
	if (seq.seq.count("N") + seq.seq.count("n") + seq.seq.count("?") + seq.seq.count("-")) / len(seq.seq) < 0.5:
		tmp=species[species['sseqid'] == seq.id]['species'].tolist()
		tmp="".join(tmp).replace(" ","_")
		seq.id=seq.id+"_"+tmp
		sequences2.append(seq)

out = open("Concatenated.phy","w")
out.write(' '+str(len(sequences2))+' '+str(len(sequences2[0].seq)))
for seq in sequences2:
	out.write("\n")	
	out.write(seq.id+" "+str(seq.seq))

out.close()

############################################### RUNNING CONCATENATED PHYLOGENY

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running Concentated Phylogeny :::\n")
sp.call("iqtree -s Concatenated.phy -bb 1000 -seed 12345", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
tree = Phylo.read("Concatenated.phy.contree", 'newick')
tree.root_at_midpoint()
Phylo.write(tree,"Concatenated.phy.contree2","newick")
Phylo.draw_ascii(tree)

sp.call('rm *.trim.*', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: FINISHED :::\n")
