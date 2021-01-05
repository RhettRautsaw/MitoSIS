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
import pandas as pd
from dfply import *
import glob

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
  ___                      _     ___ ___   ___ _           _           
 | _ ) __ _ _ _ __ ___  __| |___|_ _|   \ / __| |_  ___ __| |_____ _ _ 
 | _ \/ _` | '_/ _/ _ \/ _` / -_)| || |) | (__| ' \/ -_) _| / / -_) '_|
 |___/\__,_|_| \__\___/\__,_\___|___|___/ \___|_||_\___\__|_\_\___|_|  
                                                                       
                                                    
                     ▄▄▄▄▄▄▄   ▄ ▄▄  ▄  ▄▄ ▄▄▄▄▄▄▄  
                     █ ▄▄▄ █ ▀ ▀ ▄▀█▀▀ ▄   █ ▄▄▄ █  
                     █ ███ █ ▀█▀ ▀▀▀▀▄█▄▄█ █ ███ █  
                     █▄▄▄▄▄█ █▀▄▀█ █▀▄ █ █ █▄▄▄▄▄█  
                     ▄▄▄▄▄ ▄▄▄█▀█▄▄▄▄▄▄  ▀▄ ▄ ▄ ▄   
                      ▄▄█  ▄▄▀ █ █ ▄▄█▀▀▀▀ ▀▀█   ▀  
                      ▄▀▀ ▄▄ ▀▄▄▀▀ ██ ▄▀▀▄    █▄▀   
                     ▀██ ▄▄▄█ ▀ ▀ ▀ ▄▀ ▀▀██▀▀█▄▄ ▀  
                      ▄█  ▄▄▄█▀▀ ▄ █▄ ▄ █▀█▀█▀▄▄▀   
                     █▀  ▄▄▄▀█▀ ▀  ▄██▄▄▀▄▀███ █ ▀  
                     █ ██▀▀▄█   ██ ▄ ▀▄▀▄███▄▄ ▄█▄  
                     ▄▄▄▄▄▄▄ █ █▀ ▄ ▄ ▀█▄█ ▄ ███▀▀  
                     █ ▄▄▄ █ ▄ ▀▀█▄███▄▀▀█▄▄▄█ ▄█   
                     █ ███ █ █ ▄▀ ▀ ▄▀ ▀██▄▄▄▄██▄▀  
                     █▄▄▄▄▄█ █▀█ ▄▄▄██▄ ▄▀▀▄█▄▀▄▀   
                                                    
This script is a wrapper for MitoZ and MITGARD and is designed extract vertebrate mitochondrial genomes from raw or trimmed sequencing data. It will then blast the resulting mitogenome or barcoding genes (e.g., CYTB, COX1, ND4, 16S, etc.) to check for sample mislabeling and produce a phylogeny using IQTREE. 

:: PIPELINE ::
1. Map fastq reads to reference fasta using bwa
2. Keep reads that successfully mapped using samtools
3. Assemble mitogenome using MitoZ
4. If MitoZ fails, identify the best reference sequence
5. Use MITGARD to assemble mitogenome
	5.1 Use MitoZ to annotate mitogenome (if not already done)
6. Extract protein coding/barcoding genes
7. Blast mitogenome or genes to reference database
8. Export results and sequences
9. Align sequences and build phylogeny 

:: EXAMPLE ::
BarcodeIDChecker.py -f1 sample_F.fastq.gz -f2 sample_R.fastq.gz -r 2020-09_GenbankSnakeMito.fasta -o sample -c 16 -M 55G

:: CITE :: 
https://github.com/RhettRautsaw/BarcodeIDChecker\n\n""")

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
parser.add_argument("--version", action='version', version='BarcodeIDChecker v2')
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

#reference_name_gb="/zfs/venom/Rhett/Data/2020-10_GenbankSnakeMito/2020-10_SnakeMito.gb"
#reference_name = reference_name_gb.split(".gb")[0]+".fasta"
#references = list(SeqIO.parse(reference_name,"fasta"))

fastq1_name = os.path.abspath(args.fastq1.name)
fastq2_name = os.path.abspath(args.fastq2.name)
reference_name_gb = os.path.abspath(args.reference.name)
reference_name = reference_name_gb.split(".gb")[0]+".fasta"
output = args.output
num_threads = args.cpu
memory=args.memory

print("""
  ___                      _     ___ ___   ___ _           _           
 | _ ) __ _ _ _ __ ___  __| |___|_ _|   \ / __| |_  ___ __| |_____ _ _ 
 | _ \/ _` | '_/ _/ _ \/ _` / -_)| || |) | (__| ' \/ -_) _| / / -_) '_|
 |___/\__,_|_| \__\___/\__,_\___|___|___/ \___|_||_\___\__|_\_\___|_|  
                                                                       
                                                    
                     ▄▄▄▄▄▄▄   ▄ ▄▄  ▄  ▄▄ ▄▄▄▄▄▄▄  
                     █ ▄▄▄ █ ▀ ▀ ▄▀█▀▀ ▄   █ ▄▄▄ █  
                     █ ███ █ ▀█▀ ▀▀▀▀▄█▄▄█ █ ███ █  
                     █▄▄▄▄▄█ █▀▄▀█ █▀▄ █ █ █▄▄▄▄▄█  
                     ▄▄▄▄▄ ▄▄▄█▀█▄▄▄▄▄▄  ▀▄ ▄ ▄ ▄   
                      ▄▄█  ▄▄▀ █ █ ▄▄█▀▀▀▀ ▀▀█   ▀  
                      ▄▀▀ ▄▄ ▀▄▄▀▀ ██ ▄▀▀▄    █▄▀   
                     ▀██ ▄▄▄█ ▀ ▀ ▀ ▄▀ ▀▀██▀▀█▄▄ ▀  
                      ▄█  ▄▄▄█▀▀ ▄ █▄ ▄ █▀█▀█▀▄▄▀   
                     █▀  ▄▄▄▀█▀ ▀  ▄██▄▄▀▄▀███ █ ▀  
                     █ ██▀▀▄█   ██ ▄ ▀▄▀▄███▄▄ ▄█▄  
                     ▄▄▄▄▄▄▄ █ █▀ ▄ ▄ ▀█▄█ ▄ ███▀▀  
                     █ ▄▄▄ █ ▄ ▀▀█▄███▄▀▀█▄▄▄█ ▄█   
                     █ ███ █ █ ▄▀ ▀ ▄▀ ▀██▄▄▄▄██▄▀  
                     █▄▄▄▄▄█ █▀█ ▄▄▄██▄ ▄▀▀▄█▄▀▄▀   
                                                    
""")
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting BarcodeIDChecker...")
CWD = os.getcwd()
print("\tForward Reads -> "+fastq1_name)
print("\tReverse Reads -> "+fastq2_name)
print("\tReference Database -> "+reference_name)
print("\tOutput -> " + CWD + "/BarcodeIDChecker_results/"+output+"*")
print("\tNumber of CPU -> "+str(num_threads))
print("\tAmount of memory -> "+memory+"\n")

os.mkdir("BarcodeIDChecker_results")
os.chdir('BarcodeIDChecker_results')

############################################### PARSE GENBANK

if os.path.isfile(reference_name + ".sp"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Genbank to Fasta conversion already done :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting Genbank to Fasta :::\n")
	sp.call('rm ' + reference_name, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	species=[]
	count=1
	gb = open(reference_name_gb, 'r')
	for gb_record in SeqIO.parse(gb,'genbank'):
		acc = gb_record.id
		org = gb_record.annotations['organism']
		acc_org = [acc,org]
		species.append(acc_org)
		with open(reference_name, "a") as output_handle:
			tmp=SeqIO.write(gb_record, output_handle, "fasta")
			count = 1 + count
	
	species_df = pd.DataFrame(species)
	species_df.to_csv(reference_name+'.sp', index=False, header=False, sep="\t")
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converted %i Genbank records to Fasta :::\n" % count)

references = list(SeqIO.parse(reference_name,"fasta"))
species = pd.read_csv(reference_name+'.sp',sep="\t", names=['sseqid','species']) 

############################################### RUN BWA

if os.path.isfile(reference_name + ".amb"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: bwa index already run :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running bwa index :::\n")
	sp.call('bwa index ' + reference_name, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running bwa mem :::\n")
sp.call("bwa mem -t " + str(num_threads) + " " + reference_name + " " + fastq1_name + " " + fastq2_name + " > tmp.sam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### EXTRACT MAPPED READS WITH SAMTOOLS

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Sorting/converting sam files :::\n")
sp.call("samtools sort tmp.sam > tmp.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Pulling out reads that successfully mapped  :::\n")
sp.call("samtools view -b -F 4 tmp.bam > tmp2.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting reads to fastq format :::\n")
sp.call("samtools fastq -1 tmp_F.fq -2 tmp_R.fq -n tmp2.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### RUN MITOZ

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running MitoZ assembly :::\n")
sp.call("MitoZ.py all2 --genetic_code auto --clade Chordata --outprefix " + output + " --thread_number " + str(num_threads) + " --fastq1 tmp_F.fq --fastq2 tmp_R.fq --run_mode 2", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### IF MITOZ FAILED, RUN MITGARD

if os.path.isdir(output+".result"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: MitoZ ran successfully :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: MitoZ failed, trying MITGARD :::\n")
	sp.call("rm -rf " + output + ".tmp/", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("samtools index tmp2.bam", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	############################################### IDENTIFY BEST REFERENCE SEQUENCE BY NUMBER OF MAPPED READS
	
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Identifying best reference sequence :::\n")
	
	sp.call("samtools idxstats tmp2.bam | awk '$2 > 10000' | sort -k3nr > alternate_references.txt", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("cut -f1 alternate_references.txt | head -1 > best_reference.txt", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	
	with open('best_reference.txt', 'r') as file:
		best_reference = file.read().replace('\n', '')
	
	ref_seq = []
	for seq in references:
		if seq.id == best_reference:
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
	sp.call("MitoZ.py annotate --genetic_code auto --clade Chordata --outprefix " + output + " --thread_number " + str(num_threads) + " --fastafile tmp_mitogenome.fa", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

if os.path.isdir(output + ".result"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Moving onto BLAST :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: MitoZ annotation failed :::\n")
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Genome assembly too fragmented :::\n")
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: We will blast what MITGARD was able to assemble :::\n")

############################################### CREATE BLAST DATABASE FROM REFERENCE SEQUENCES

if os.path.isfile(reference_name + ".nin"):
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: makeblastdb already run :::\n")
else:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running makeblastdb :::\n")
	sp.call('makeblastdb -in ' + reference_name + ' -dbtype nucl', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

############################################### BLAST ANNOTATED REGIONS OR FULL MITOGENOME (IF ANNOTATION FAILS) TO REFERENCE DATABASE

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running blast :::\n")
if os.path.isdir(output + ".result"):
	sp.call('mv '+output+'.result mitoz.result', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('cat mitoz.result/' + output + '.cds mitoz.result/' + output + '.rrna > blast_query.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("perl -pi -e 's/>.*?;/>"+output+";/g' blast_query.fasta", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('cp mitoz.result/*most_related_species.txt mitoz_most_related_species.txt', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('cp mitoz.result/'+output+'.fasta '+output+'_mitogenome.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("perl -pi -e 's/>.*?;/>"+output+";/g' "+output+"_mitogenome.fasta", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('blastn -query blast_query.fasta -db ' + reference_name + ' -outfmt "6 qseqid sseqid stitle pident evalue bitscore sseq" -num_threads ' + str(num_threads) + ' -max_target_seqs 50 -qcov_hsp_perc 90 -evalue 0.0001 -out blast_results.tab', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	query = list(SeqIO.parse("blast_query.fasta","fasta"))
else:
	sp.call("sed -i 's/>scaffold1/>"+output+"/g' tmp_mitogenome.fa", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('mv tmp_mitogenome.fa '+output+'_mitogenome.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call('blastn -query '+output+'_mitogenome.fasta -db ' + reference_name + ' -outfmt "6 qseqid sseqid stitle pident evalue bitscore sseq" -num_threads ' + str(num_threads) + ' -max_target_seqs 50 -evalue 0.0001 -out blast_results.tab', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	query = list(SeqIO.parse(output+'_mitogenome.fasta',"fasta"))

############################################### EXTRACT BLAST MATCH SEQUENCES

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Extracting BLAST results for Phylogenetics :::\n")
header_list = ["qseqid","sseqid","stitle","pident","evalue","bitscore","sseq"]
results = pd.read_csv("blast_results.tab",sep='\t',names=header_list)
results = pd.merge(results, species, on ='sseqid', how ='left')

results2=results >> group_by(X.species) >> summarize(mean_pident = X.pident.mean())
results2=results2.sort_values("mean_pident", ascending=False)
tfile = open('blast_results.summary', 'w')
tfile.write(results2.to_string(index=False))
tfile.close()
#results2.to_csv("blast_results_summary.csv",index=False)

sp.call("sort -t$'\t' -k4 -nr blast_results.tab > tmp.tab", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call("mv tmp.tab blast_results.tab", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
sp.call('rm -rf tmp* *.tmp mapped bowtie_index align.sam consensus.mfa.fasta', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

os.mkdir("Phylogenetics")
os.chdir("Phylogenetics")

if os.path.isdir("../mitoz.result"):
	files=[]
	for i in results.qseqid.unique():
		gene=i.split(";")[1]
		files.append(gene+'.fasta')
		print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Extracting " + gene + " BLAST matches :::\n")
		
		############################################### CONVERTING PANDAS DATAFRAME (BLAST OUTFMT 6) TO FASTA
		
		ref_seqs=results >> mask(X.qseqid==i) >> select(X.stitle, X.sseq) >> mutate(stitle=">"+X.stitle)
		fna=[]
		for index in ref_seqs.index:
			fna.append(ref_seqs['stitle'][index].replace(" ","_"))
			fna.append(ref_seqs['sseq'][index].replace("-",""))
		
		with open(gene+'.fasta', 'w') as f:
			for item in fna:
				bytes=f.write('%s\n' % item)
		
		############################################### COMBINING NEW FASTA WITH QUERY SEQUENCE
		
		ref_seqs = list(SeqIO.parse(gene+'.fasta',"fasta"))
		
		gene_db = []
		for seq in query:
			if seq.id == i:
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

for i in files:
	print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Aligning, Trimming, and Inferring Phylogeny for " + i + " :::\n")
	sp.call("sed -i 's/ /_/g' "+i, shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("mafft "+i+" > "+i+".aln", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("trimal -in "+i+".aln -out "+i+".trim -automated1", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("iqtree -s "+i+".trim -bb 1000 -seed 12345", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	if os.path.isfile(i+".trim.contree"):
		tree = Phylo.read(i+".trim.contree", 'newick')
		tree.root_at_midpoint()
		Phylo.write(tree,i+".contree","newick")
		Phylo.draw_ascii(tree)
		sp.call('mv '+i+'.trim.iqtree '+i+'.iqtree',shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

sp.call('rm *.trim.*', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: FINISHED :::\n")
