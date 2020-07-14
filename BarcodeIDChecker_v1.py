#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='Barcode ID Checker\n\n\
This script is a wrapper for many other programs and is designed extract the mitochondrial genome from raw or trimmed sequencing data and blast the resulting protein coding sequences or barcoding genes (e.g., CYTB, COX1, ND4, etc.) to check for sample mislabeling. \n\
The second part of this script will use the assembled mitogenomes to infer a phylogeny using HOMBLOCKS and IQTREE. \n\n\
:: INPUT :: \n\
--assemble \n\
-f1, f2 \t Fastq files containing the raw/trimmed reads \n\
-r \t\t A fasta file containing reference sequences to map reads to prior to assembly and blast to following assembly. \n\n\
--phylo \n\
-s \t\t A txt list containing the names of all samples on which assemble has been run. \n\n\
:: PIPELINE ::\n\
--assemble \n\
- Map fastq reads to reference fasta using bwa \n\
- Keep reads that successfully mapped using samtools \n\
- Assemble mitogenome using MitoZ \n\
- Extract protein coding/barcoding genes \n\
- Blast genes to reference fasta and determine top blast hit for taxa ID \n\n\
--phylo \n\
- Align mitogenomes using HomBlocks \n\
- Infer phylogeny using IQTREE \n\n\
:: EXAMPLES :: \n\
BarcodeIDChecker_v1.py --assemble -f1 {}_F.fq -f2 {}_R.fq -r mito_ref.fasta -o {}\n\
BarcodeIDChecker_v1.py --phylo -s samples.txt \n\n\
::CITE:: \n\
https://github.com/reptilerhett/Bioinformatic-Scripts\n\n')

###############################################

parser.add_argument("--assemble",
					action="store_true",
					help="Perform mitogenome assembly and species ID check.")
parser.add_argument("--phylo",
					action="store_true",
					help="Infer mitogenome phylogenetic tree. \nWARNING: --phylo can only be run after --assemble is run")

parser.add_argument("-f1","--fastq1",
					type=argparse.FileType('r+'),
					help="REQUIRED for --assemble: Fastq read pair 1 (forward)")
parser.add_argument("-f2","--fastq2",
					type=argparse.FileType('r+'),
					help="REQUIRED for --assemble: Fastq read pair 2 (reverse)")
parser.add_argument("-r","--reference",
					type=argparse.FileType('r+'),
					help="REQUIRED for --assemble: Fasta reference database (I recommend downloading all mitochondrial data for your clade of interest from GenBank, e.g. snakes)")
parser.add_argument("-o","--output",
					type=str,
					default='ZZZ',
					help="OPTIONAL for --assemble: Prefix for output files. Default is 'ZZZ'")

parser.add_argument("-s","--samples",
					type=argparse.FileType('r+'),
					help="REQUIRED for --phylo: List of samples that have had Step 1 run on them.")

parser.add_argument("-t","--num_threads",
					type=int,
					default=8,
					help="OPTIONAL for both --assemble and --phylo: Number of threads. Default is 8")

parser.add_argument("--version", action='version', version='BarcodeIDChecker_v3')
args=parser.parse_args()

###############################################

#fastq1_name="/zfs/venom/ajmason/Active/MS_Phylogenomics_of_Bothrops_2019/data/02_trimming/I0771_Agkistrodon_conanti/I0771_Agkistrodon_conanti_F_val_1.fq.gz"
#fastq2_name="/zfs/venom/ajmason/Active/MS_Phylogenomics_of_Bothrops_2019/data/02_trimming/I0771_Agkistrodon_conanti/I0771_Agkistrodon_conanti_R_val_2.fq.gz"
#reference_name="../mito_ref.fasta"
#output="I0771_Agkistrodon_conanti"
#num_threads=8

output = args.output
num_threads = args.num_threads

###############################################

if args.assemble:
	print("\n::: AVENGERS ASSEMBLE! :::\n")

	os.mkdir("BarcodeIDChecker_results")

	fastq1_name = args.fastq1.name
	fastq2_name = args.fastq2.name
	reference_name = args.reference.name
	references = list(SeqIO.parse(reference_name,"fasta"))

	if os.path.isfile(reference_name + ".amb"):
		print("\n::: bwa index already run :::\n")
	else:
		print("\n::: Running bwa index :::\n")
		command = 'bwa index ' + reference_name
		subprocess.call(command, shell=True)

	print("\n::: Running bwa mem :::\n")
	command = "bwa mem -t " + str(num_threads) + " " + reference_name + " " + fastq1_name + " " + fastq2_name + " > BarcodeIDChecker_results/tmp.sam"
	subprocess.call(command, shell=True)

	os.chdir('BarcodeIDChecker_results')

	print("\n::: Sorting/converting sam files :::\n")
	command = "samtools sort tmp.sam > tmp.bam"
	subprocess.call(command, shell=True)

	print("\n::: Pulling out reads that successfully mapped  :::\n")
	command = "samtools view -b -F 4 tmp.bam > tmp2.bam"
	subprocess.call(command, shell=True)

	print("\n::: Converting reads to fastq format :::\n")
	command = "samtools fastq -1 " + output + "_F.fq -2 " + output + "_R.fq -n tmp2.bam"
	subprocess.call(command, shell=True)

	print("\n::: Running MitoZ assembly :::\n")
	command = "MitoZ.py all2 --genetic_code auto --clade Chordata --outprefix " + output + " --thread_number " + str(num_threads) + " --fastq1 " + output + "_F.fq --fastq2 " + output + "_R.fq --run_mode 2"
	subprocess.call(command, shell=True)

	if os.path.isfile("./*.result"):
		print("\n::: MitoZ ran successfully :::\n")
	else:
		print("\n::: MitoZ did not run successfully :::\n")
		print("::: Using MITGARD instead :::\n")

		subprocess.call("rm -rf " + output + ".tmp/", shell=True)

		subprocess.call("samtools index tmp2.bam", shell=True)

		print("\n::: Identifying best reference sequence :::\n")

		command = "samtools idxstats tmp2.bam | awk '$2 > 10000' | sort -k3nr | cut -f1 | head -1 > best_reference.txt"
		subprocess.call(command,shell=True)

		command = "samtools idxstats tmp2.bam | awk '$2 > 10000' | sort -k3nr > alternate_references.txt"
		subprocess.call(command,shell=True)

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

		subprocess.call("perl -pi -e 's/ .*//g' best_reference.fasta", shell=True)

		print("\n::: Running MITGARD :::\n")
		command = "MITGARD.py -s tmp -1 " + fastq1_name + " -2 " + fastq2_name + " -R best_reference.fasta -c " + str(num_threads)
		subprocess.call(command,shell=True)

		subprocess.call("perl -pi -e 's/>.*/>scaffold1/g' tmp_mitogenome.fa", shell=True)

		print("\n::: Annotating MITGARD mitogenome with MitoZ :::\n")
		command = "MitoZ.py annotate --genetic_code auto --clade Chordata --outprefix " + output + " --thread_number " + str(num_threads) + " --fastafile tmp_mitogenome.fa"
		subprocess.call(command,shell=True)

	if os.path.isfile("../" + reference_name + ".nin"):
		print("\n::: makeblastdb already run :::\n")
	else:
		print("\n::: Running makeblastdb :::\n")
		command = 'makeblastdb -in ../' + reference_name + ' -dbtype nucl'
		subprocess.call(command, shell=True)

	print("\n::: Running blast :::\n")
	if os.path.isdir(output + ".result"):
		subprocess.call('cat ' + output + '.result/' + output '.cds ' + output + '.result/' + output + '.rrna > blast_query.fasta', shell=True)
		
		command = 'blastn -query blast_query.fasta -db ../' + reference_name + ' -outfmt "6 qseqid stitle pident evalue bitscore" -num_threads ' + str(num_threads) + ' -max_target_seqs 10 -evalue 0.0001 -out ' + output + "_blast.tab"
		subprocess.call(command, shell=True)
	else:
		print("\n::: MitoZ annotation failed :::\n")
		print("\n::: Genome assembly too fragmented :::\n")
		print("\n::: We will blast what MITGARD was able to assemble :::\n")
		command = 'blastn -query tmp_mitogenome.fa -db ../' + reference_name + ' -outfmt "6 qseqid stitle pident evalue bitscore" -num_threads ' + str(num_threads) + ' -max_target_seqs 10 -evalue 0.0001 -out ' + output + "_blast.tab"
		subprocess.call(command, shell=True)

		os.mkdir(output + ".result")
		subprocess.call('mv tmp_mitogenome.fa ' + output + '.result/' + output + '_mitogenome.fasta', shell=True)

	subprocess.call("sort -t$'\t' -k3 -nr " + output + "_blast.tab > tmp.tab", shell=True)

	subprocess.call("mv tmp.tab " + output + "_blast.tab", shell=True)

	command = 'rm -rf tmp* *.tmp mapped bowtie_index align.sam consensus.mfa.fasta'
	subprocess.call(command, shell=True)

if args.phylo:
	print("\n::: AVENGERS INFER PHYLOGENY! :::\n")

	samples_name = args.samples.name

#	os.mkdir("HOMBLOCKS")

#	command = "perl HomBlocks.pl --align --path=/public/home/mgb217/HomBlocks/Xenarthrans/fasta/ -out_seq=Xenarthrans.output.fasta  --mauve-out=Xenarthrans.mauve.out


print("\n::: FINISHED :::\n")
