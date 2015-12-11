#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net

import os
import sys
import re
import glob
import shutil
import argparse
import ConfigParser
import commands
import subprocess
from Bio import SeqIO

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Input %%%


# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Map raw reads against reference sequence and phase the reads into two separate alleles. Then produce consensus sequence for each allele.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--reads',
		required=True,
		action=CompletePath,
		default=None,
		help='Call the folder that contains the trimmed reads, organized in a separate subfolder for each sample. The name of the subfolder has to start with the sample name, delimited with an underscore [_].'
	)	
	parser.add_argument(
		'--alignments',
		required=True,
		action=CompletePath,
		default=None,
		help='The folder containing exon-/gene-alignments of your target loci.'
	)
	parser.add_argument(
		'--config',
		required=True,
		help='A configuration file containing the full paths to the following programs: samtools, clc-assembly-cell, bcftools, vcfutils, picard'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed.'
	)
	parser.add_argument(
		'--no_duplicates',
		action='store_true',
		default=False,
		help='Use this flag if you want to clean the mapped reads from all duplicates with Picard.'
	)
	parser.add_argument(
		'--conservative',
		action='store_true',
		default=False,
		help='Use this flag if you want to discard all base calls with limited certainty (covered by <3 reads). This will produce the ambiguity character "N" instead of that potential base call in the final sequence.'
	)	
	parser.add_argument(
		'--cores',
		type=int,
		default=1,
		help='For parallel processing you can choose the number of cores you want CLC to run on.'
	)
	return parser.parse_args()

# Preparation for calling input variables and files
args = get_args()
conf = ConfigParser.ConfigParser()
conf.optionxform = str
conf.read(args.config)

# Collect all program paths from control file
paths = conf.items('paths')

# Find each program name in list and define variable
samtools = ""	
for i in paths:
		if "samtools" in i:
			samtools = i[1]	
bcftools = ""	
for i in paths:
		if "bcftools" in i:
			bcftools = i[1]	
vcfutils = ""	
for i in paths:
		if "vcfutils" in i:
			vcfutils = i[1]
picard = ""	
for i in paths:
		if "picard" in i:
			picard = i[1]
clc = ""	
for i in paths:
		if "clc" in i:
			clc = i[1]
emboss = ""
for i in paths:
		if "emboss" in i:
			emboss = i[1]
	
# Specify the different program commands
mark_duplicates = os.path.join(picard,"MarkDuplicates.jar")
cons = os.path.join(emboss,"bin/cons")
clc_mapper = os.path.join(clc,"clc_mapper")
clc_cas_to_sam = os.path.join(clc,"clc_cas_to_sam")
clc_assembler = "%s/clc_assembler" %clc

# Set working directory
out_dir = args.output
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
	
# Get other input variables
alignments = args.alignments
reads = args.reads


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Functions %%%


def create_reference_fasta(sample_id):
	print "Creating reference library for %s .........." %sample_id
#	get the sequence header with the correct fasta id and extract sequence
#	store these sequences in separate fasta file for each locus at out_dir/reference_seqs/sample_id
#	header of sequence remains the locus name
#	remove all "-" and "?" in sequence (not in header!!)
	sample_reference_folder = os.path.join(reference_folder,sample_id)
	if not os.path.exists(sample_reference_folder):
		os.makedirs(sample_reference_folder)
	file_list = [fn for fn in os.listdir(alignments) if fn.endswith(".fasta")]
	for fasta_alignment in file_list:
		locus_id = fasta_alignment.replace(".fasta", "")
		sample_reference_fasta = os.path.join(sample_reference_folder,fasta_alignment)
		fasta_sequences = SeqIO.parse(open("%s/%s" %(alignments,fasta_alignment)),'fasta')
		outfile = open(sample_reference_fasta, 'w')
		for fasta in fasta_sequences:
			if fasta.id == sample_id:
				outfile.write(">%s\n%s\n" %(locus_id,fasta.seq))
		outfile.close()
		remove_chars = "sed -i 's/[-?]//g' %s" %sample_reference_fasta
		os.system(remove_chars)
	reference = os.path.join(sample_reference_folder,"joined_fasta_library.fasta")
	join_fastas = "cat %s/*.fasta > %s" %(sample_reference_folder,reference)
	os.system(join_fastas)
	return reference


	
	


def mapping(forward,backward,reference,sample_id,sample_output_folder):
	print "Mapping.........."
	cas = "%s/%s.cas" %(sample_output_folder,sample_id)
	command1 = "%s -o %s -d %s -q -p fb ss 100 1000 -i %s %s -l 0.9 -s 0.9 --cpus %d" %(clc_mapper,cas,reference,forward,backward,args.cores)
	os.system(command1)
	
	print "Converting to bam.........."
	bam = "%s/%s.bam" %(sample_output_folder,sample_id)
	command2 = "%s -a %s -o %s -f 33 -u" %(clc_cas_to_sam,cas,bam)
	os.system(command2)

	print "Sorting bam.........."
	sorted = "%s/%s.sorted" %(sample_output_folder,sample_id)
	command3 = "%s sort %s %s" %(samtools,bam,sorted)
	os.system(command3)

	print "Indexing bam.........."	
	sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
	command4 = "%s index %s" %(samtools,sorted_bam)
	os.system(command4)
	
	print "Removing obsolete files.........."
	command5 = "rm %s %s" %(cas,bam)
	os.system(command5)
	
	return sorted_bam
	

def clean_with_picard(sample_output_folder,sample_id,sorted_bam):
	picard_folder = "%s/picard" %sample_output_folder
	if not os.path.exists(picard_folder):
		os.makedirs(picard_folder)
	picard_log_folder = "%s/log" %picard_folder
	if not os.path.exists(picard_log_folder):
		os.makedirs(picard_log_folder)	
	picard_out = "%s/%s_no_dupls_sorted.bam" %(picard_folder,sample_id)
	dupl_log = "%s/%s_dupls.log" %(picard_log_folder,sample_id)
	run_picard = [
		"java",
		"-jar",
		mark_duplicates,
		"I=%s" %sorted_bam,
		"O=%s" %picard_out,
		"M=%s" %dupl_log,
		"REMOVE_DUPLICATES=true", 
		"VALIDATION_STRINGENCY=LENIENT"
	]
	try:
		print "Removing duplicate reads with Picard.........."		
		with open(os.path.join(picard_log_folder, "picard_screen_out.txt"), 'w') as log_err_file:		
			pi = subprocess.Popen(run_picard, stderr=log_err_file)
			pi.communicate()
		print "Duplicates successfully removed."
	except:
		print "Running Picard caused an error. Please check your picard path specification in the control file."
		sys.exit()
		
	print "Indexing Picard-cleaned bam.........."	
	index_picard_bam = "%s index %s" %(samtools,picard_out)
	os.system(index_picard_bam)
	return picard_out
	
	
def phase_bam(sorted_bam_file,sample_output_folder,reference):
	# Phasing:
	bam_basename = re.sub('.bam$', '', sorted_bam_file)
	split_sample_path = re.split("/",sorted_bam_file)
	split_file_name = split_sample_path[-1]
	phasing_file_base_pre = re.sub('.sorted.bam$', '', split_file_name)
	
	phasing_out_dir = "%s/phased" %(sample_output_folder)
	if not os.path.exists(phasing_out_dir):
		os.makedirs(phasing_out_dir)
	phasing_basename = "%s/%s_allele" %(phasing_out_dir,phasing_file_base_pre)
	phasing_cmd = [
		samtools,
		"phase", 
		"-A",
		"-F",
		"-Q",
		"20",
		"-b",
		phasing_basename,
		sorted_bam_file
	]
	try:
		print "Phasing bam file.........."
		with open(os.path.join(phasing_out_dir, "phasing_screen_out.txt"), 'w') as phasing_screen:		
			ph = subprocess.Popen(phasing_cmd, stdout=phasing_screen)
			ph.communicate()
		print "Phasing completed."
	except:
		print "Phasing unsuccessful. Script terminated."
		sys.exit()
	
	allele_0_file = "%s.0.bam" %phasing_basename
	allele_1_file = "%s.1.bam" %phasing_basename
	allele_0_sorted_base = "%s/%s_sorted_allele_0" %(phasing_out_dir,phasing_file_base_pre)
	allele_1_sorted_base = "%s/%s_sorted_allele_1" %(phasing_out_dir,phasing_file_base_pre)
	allele_0_sorted_file = "%s.bam" %allele_0_sorted_base
	allele_1_sorted_file = "%s.bam" %allele_1_sorted_base
	
	# Sorting phased bam files:
	sort_phased_0 = "%s sort %s %s" %(samtools,allele_0_file,allele_0_sorted_base)
	sort_phased_1 = "%s sort %s %s" %(samtools,allele_1_file,allele_1_sorted_base)
	os.system(sort_phased_0)
	os.system(sort_phased_1)
	
	# Creating consensus sequences from bam-files
	print "Creating consensus sequences from bam-files.........."
	make_cons_0 = "%s mpileup -u -f %s %s | %s view -cg - | %s vcf2fq > %s_0.fq" %(samtools,reference,allele_0_sorted_file,bcftools,vcfutils,phasing_basename)
	make_cons_1 = "%s mpileup -u -f %s %s | %s view -cg - | %s vcf2fq > %s_1.fq" %(samtools,reference,allele_1_sorted_file,bcftools,vcfutils,phasing_basename)
	make_cons_unphased = "%s mpileup -u -f %s %s | %s view -cg - | %s vcf2fq > %s.fq" %(samtools,reference,sorted_bam_file,bcftools,vcfutils,bam_basename)
	os.system(make_cons_0)
	os.system(make_cons_1)
	os.system(make_cons_unphased)
	
	# Converting fq into fasta files
	make_fasta_cmd_0 = "seqtk seq -a %s_0.fq > %s_0.fasta" %(phasing_basename,phasing_basename)
	make_fasta_cmd_1 = "seqtk seq -a %s_1.fq > %s_1.fasta" %(phasing_basename,phasing_basename)
	make_fasta_cmd_unphased = "seqtk seq -a %s.fq > %s.fasta" %(bam_basename,bam_basename)
	os.system(make_fasta_cmd_0)
	os.system(make_fasta_cmd_1)
	os.system(make_fasta_cmd_unphased)
	
	# Cleaning up output directory 
	output_files = [val for sublist in [[os.path.join(i[0], j) for j in i[2]] for i in os.walk(sample_output_folder)] for val in sublist]
	fasta_dir = "%s/final_fasta_files" %sample_output_folder
	if not os.path.exists(fasta_dir):
		os.makedirs(fasta_dir)
	# check the names to make sure we're not deleting something improperly
	try:
		assert "%s.fasta" %bam_basename in output_files
	except:
		raise IOError("Output-files were not created properly.")
	for file in output_files:
		if file in ("%s.fasta" %bam_basename, "%s_0.fasta" %phasing_basename, "%s_1.fasta" %phasing_basename):
			shutil.move(file,fasta_dir)

	if args.conservative:
		# Clean up the final fasta alignments and replace all uncertain base-calls (non-capitalized letters) with "N"
		replace_uncertain_base_calls = "for fasta in $(ls %s/*.fasta); do sed -i -e '/>/! s=[actgn]=N=g' $fasta; done" %fasta_dir		
		os.system(replace_uncertain_base_calls)
	
	# Standardize all the different ambiguity code-bases with N
	standardize_all_ambiguities = "for fasta in $(ls %s/*.fasta); do sed -i -e '/>/! s=[ywrksmYWRKSM]=N=g' $fasta; done" %fasta_dir
	os.system(standardize_all_ambiguities)
	# Remove the unnecessary .fq files
	remove_fq_file = "rm %s/*.fq" %sample_output_folder
	remove_fq_file_2 = "rm %s/*/*.fq" %sample_output_folder
	if args.no_duplicates:
		os.system(remove_fq_file_2)
	else: 
		os.system(remove_fq_file)
	return fasta_dir
		
		
def edit_fasta_headers(allele_fastas,sample_id):
	cmd0 = "sed -i -e 's/>\(.*\)/&_%s_0 |&_phased/g' %s/*allele_0.fasta" %(sample_id,allele_fastas)
	os.system(cmd0)
	cmd1 = "sed -i -e 's/>\(.*\)/&_%s_1 |&_phased/g' %s/*allele_1.fasta" %(sample_id,allele_fastas)
	os.system(cmd1)
	cmd_unphased = "sed -i -e 's/>\(.*\)/&_%s_hom |&/g' %s/*sorted.fasta" %(sample_id,allele_fastas)
	os.system(cmd_unphased)	
	cmd_final = "for allele in $(ls %s/*.fasta); do sed -i 's/|>/|/g' $allele; done" %allele_fastas
	os.system(cmd_final)
	
	
def join_allele_fastas():
	final_merging = "for folder in $(find %s -type d -name '*_remapped'); do cat $folder/final_fasta_files/*allele*; done > %s/joined_allele_sequences_all_samples.fasta" %(out_dir,out_dir)
	os.system(final_merging)
	capitalize = "sed -i -e '/>/! s=[actgn]=\U&=g' %s/joined_allele_sequences_all_samples.fasta" %out_dir
	os.system(capitalize)
	

def manage_homzygous_samples(fasta_dir, sample_id):	
	fasta_sequences = SeqIO.parse(open("%s/%s.sorted.fasta" %(fasta_dir,sample_id)),'fasta')
	with open('%s/%s_joined_homozygous_alleles.fasta'%(fasta_dir,sample_id), 'w') as outfile:
		for fasta in fasta_sequences:
			name = re.split(" ", fasta.description)
			name[0] += "_0"
			fasta.description = " ".join(name)
			fasta.id += "_0"
			SeqIO.write(fasta, outfile, "fasta")
			name = re.split(" ", fasta.description)
			allele_1_name = re.sub("_0$", "_1", name[0])		
			name[0] = allele_1_name
			fasta.description = " ".join(name)
			allele_1_id = re.sub("_0$", "_1", str(fasta.id))	
			fasta.id = allele_1_id
			SeqIO.write(fasta, outfile, "fasta")
	outfile.close
	
	
def assembly_clc(forward,backward):
	print "De-novo assembly with CLC.........."	
	output_fasta = "%s/%s_contigs.fasta" %(sample_output_folder,sample_id)
	distance_file = "%s/%s-estimated-distances.txt" %(sample_output_folder,sample_id)
	command = "%s -o %s -q -i %s %s -e %s -p fb ss 150 800 --cpus %d" %(clc_assembler,output_fasta,forward,backward,distance_file,args.cores)
	os.system(command)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Workflow %%%

reference_folder = "%s/reference_seqs" %out_dir
if not os.path.exists(reference_folder):
	os.makedirs(reference_folder)
	


for subfolder, dirs, files in os.walk(reads):
	subfolder_path_elements = re.split("%s/" %reads, subfolder)
	if subfolder_path_elements[-1] != reads:
		sample_folder = subfolder_path_elements[-1]
		sample_id = re.sub("_clean","",sample_folder)
		reference = create_reference_fasta(sample_id)
		# Loop through each sample-folder and find read-files
		sample_output_folder = "%s/%s_remapped" %(out_dir,sample_id)
		if not os.path.exists(sample_output_folder):
			os.makedirs(sample_output_folder)
		for misc1, misc2, fastq in os.walk(subfolder):
			forward = ""
			backward = ""			
			for element in fastq:
				if sample_id in element and element.endswith("READ1.fastq"):
					forward = "%s/%s" %(subfolder,element)
				if sample_id in element and element.endswith("READ2.fastq"):
					backward = "%s/%s" %(subfolder,element)
			if forward != "" and backward != "":
				print "\n", "#" * 50
				print "Processing sample", sample_id, "\n"
				sorted_bam = mapping(forward,backward,reference,sample_id,sample_output_folder)
				if args.no_duplicates:
					sorted_bam = clean_with_picard(sample_output_folder,sample_id,sorted_bam)
				allele_fastas = phase_bam(sorted_bam,sample_output_folder,reference)

## THIS iS STILL NOT WORKING PROPERLY WHEN ONLY SINGLE FILE PRESENT:				
				
				# The following is for the case that no phased bam files were created, i.e. the individual is homozygous for all loci (happens when only looking at one locus or a very few)
				allele0 = ""
				allele1 = ""
				# testing if phasing files were created
				for file in os.listdir(allele_fastas):
					if file.endswith(".fasta"):
						if "allele_0" in file:
							allele0 = file
						if "allele_1" in file:
							allele1 = file		
				if allele0 == 0:
					manage_homzygous_samples(allele_fastas,sample_id)
					os.remove(os.path.join(allele_fastas,allele0))
					os.remove(os.path.join(allele_fastas,allele1))
				# Give fasta headers the correct format for phasing script
				edit_fasta_headers(allele_fastas,sample_id)
				#os.remove(os.path.join(allele_fastas,"%s.sorted.fasta" %sample_id))
				print "\n", "#" * 50
join_allele_fastas()

					
#			else:
#				print "\nError: Read-files for sample %s could not be found. Please check if subfolders/sample-folders are named in this pattern: 'sampleID_clean' and if the cleaned fastq files in the sample-folder end with 'READ1.fastq' and 'READ2.fastq' respectively." %sample_id
#				raise SystemExit
#	else:
#		print "\nError: Check your folder structure. The folder given at the --reads flag must contain a separate subfolder for each sample. Please check if these sample specific subfolders are named in this pattern: 'sampleID_clean' and if the cleaned fastq files in the sample-folder end with 'READ1.fastq' and 'READ2.fastq' respectively. " 
#		raise SystemExit
