#!/usr/bin/env python3
# wrapper for comparing triploidy and diploidy with microhap panels
# for the Eagle Fish Genetics Laboratory

import re
import sys
import os
import glob
import multiprocessing as mp
import subprocess
from scipy.stats import binom_test
import matplotlib.pyplot as plt


# calls short read aligner
#  executed here to allow parallel processing over samples
def alignReads(fastq, alignRef):
	stdErr = ""
	if os.path.exists(fastq):
		outfile = re.sub("\.fastq$", ".bam", fastq)
		shellOut = subprocess.run("bowtie2 --end-to-end --no-unal -N 1 --rdg 0,5 --rfg 0,5 --score-min L,0,-.76 -x " + 
			alignRef + " -U " + fastq + " | samtools view -F 276 -u | samtools sort -o " + outfile, 
			shell = True, stderr = subprocess.PIPE, universal_newlines = True)
		stdErr = shellOut.stderr
	return [fastq, stdErr]
	
# choose one poten het biallelic SNP within an indiv
#  from a microhap locus and return the read counts
def chooseSNP(counts): # counts is[["AAA", 13], ["Allele2", count], [], ...]
	if len(counts) < 2: # only one allele, so doesn't matter what you return (FOR THIS trip/dip ALGORITHM)
		return [0,0]
	else:
		al = len(counts[0][0]) # get length of allele
		allDicts = []
		for i in range(0, al, 1): # initialize empty dicts, one for each SNP in the locus
			allDicts += [{}]
		for c in counts:
			for i in range(0, al, 1):
				allDicts[i][c[0][i]] = allDicts[i].get(c[0][i], 0) + c[1] # add counts for that SNP
		# now choose a biallelic SNP, preferentially heterozygous and close to start
		for i in range(0, al, 1):
			if len(allDicts[i]) == 2:
				c = [allDicts[i][x] for x in allDicts[i]]
				if binom_test(max(c), c[0] + c[1], ((2/3)*(0.99)) + ((1/3)*.01), alternative = "greater") > 0.05: # quick screen for potential hets
					return c
		return [0,0] # none appear to be het

def Main():
	# defaults
	panel = None # -panel, panel name
	threads = len(os.sched_getaffinity(0)) # -t, number of threads to use throughout, default is max available
	bc = None # -bc, bcSplitFile.csv
	minSuccess = 0.9 # -s, minimum proportion of genotypes called to be "successful"
	prefix = "GTsEagle" # -pre, prefix for output files
	
	# get command line inputs
	flag = 1
	while flag < len(sys.argv):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == "-t":
			flag += 1
			threads = int(sys.argv[flag])
		elif sys.argv[flag] == "-panel":
			flag += 1
			panel = sys.argv[flag]
		elif sys.argv[flag] == "-bc":
			flag += 1
			bc = sys.argv[flag]
		elif sys.argv[flag] == "-s":
			flag += 1
			minSuccess = float(sys.argv[flag])
		elif sys.argv[flag] == "-pre":
			flag += 1
			prefix = sys.argv[flag]
		else:
			print("Error: option", sys.argv[flag], "not recognized.")
			return
		flag += 1
	# end while loop for command line arguments
	
	# key is panel name, value is [pos, ref fasta, ref bt2 prefix, pres abs input]
	panelDict = {
		"Sna302" : ["/home/efglserv/software/GTseq/microhapInputs/posFiles/Sna302_pos_5.txt", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/Sna302/lt_1.fasta", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/Sna302/lt_1", 
			"/home/efglserv/software/GTseq/microhapInputs/presAbsInputs/Sna302_pres_abs.txt",
			"30", "15", ".05", "-5", "100"], # min reads, min loci, binom_p_value, max LLR diploid, min LLR triploid - h assumed 1 and eps assumed .01
			
		"Svi151" : ["/home/efglserv/software/GTseq/microhapInputs/posFiles/Svi151_pos_final.txt", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/Svi151/wall_1.fasta", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/Svi151/wall_1", 
			"", # no pres abs markers for Svi151
			"30", "15", ".05", "-5", "200"],
		
		"One361" : ["/home/efglserv/software/GTseq/microhapInputs/posFiles/One361_pos_5.txt", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/One361/OneAmpRef_5.fa", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/One361/One361", 
			"", # no pres abs markers for One361
			"30", "15", ".05", "-5", "100"]
	}
		
	if panel not in panelDict:
		print("Valid panel names are:")
		for p in panelDict:
			print(p)
		return
		
	panelInfo = panelDict[panel]
	
	# align
	fastqFiles = glob.glob("*.fastq")
	bamFiles = [re.sub("\.fastq$", ".bam", x) for x in fastqFiles]

	print("Aligning reads")
	process_pool = mp.Pool(processes=threads)
	alignOut = [process_pool.apply_async(alignReads, args = [x, panelInfo[2]]) for x in fastqFiles]
	
	# output stderr (alignment stats)
	with open("alignment_output.txt", "w") as outfile:
		for i in alignOut:
			temp = i.get()
			outfile.write(temp[0] + "\n" + temp[1])
	
	process_pool.close()
	process_pool.join()

	# get allele counts with microTyper
	print("Counting reads for each allele")
	ac = subprocess.run(["mtype2", "-f"] + bamFiles + ["-p", panelInfo[0], "-r", panelInfo[1],
		"-o", "-", "-t", str(threads), "-b", "400", "--justCount"], 
		stderr = subprocess.PIPE, stdout = subprocess.PIPE, universal_newlines = True)
	ac = ac.stdout.splitlines() # Indiv, Locus, Allele, Count, removing newline characters
	ac = ac[1:] # remove header
	
	# process allele counts to choose biallelic SNPs at each locus/individual - this is really slow in R so doing it here
	print("choosing SNPs")
	with open("snp_counts.txt", "w") as outfile:
		outfile.write("\t".join(["Indiv", "Locus", "Count1", "Count2"]) + "\n")
		loc = ac[0].split("\t")[1]
		ind = ac[0].split("\t")[0]
		counts = []
		for row in ac:
			line = row.split("\t")
			if line[2] == "noReads": # this will just be output as [0,0] b/c only one "allele"
				line[3] = 0
			else:
				line[3] = int(line[3])
			if line[1] == loc and line[0] == ind:
				counts += [line[2:4]]
			else:
				counts = chooseSNP(counts)
				outfile.write("\t".join([re.sub("\.bam$", "", ind), loc, str(counts[0]), str(counts[1])]) + "\n")
				loc = line[1]
				ind = line[0]
				counts = [line[2:4]]
		# process last locus
		counts = chooseSNP(counts)
		outfile.write("\t".join([re.sub("\.bam$", "", line[0]), loc, str(counts[0]), str(counts[1])]) + "\n")
	
	# now call Rscript for tripsAndDip
	print("Inferring ploidy")
	subprocess.run(["mhTripDip.R"] + panelInfo[4:9] + [prefix])
	
	# read in and process diploids/unknown and triploids separately
	dipUnk = []
	trip = []
	with open(prefix + "_genPloidy.txt", "r") as infile:
		line = infile.readline() # skip header
		line = infile.readline() # fullName sample_name genPloidy LLR loci_used
		while line:
			sep = line.split(",")
			if sep[2] == "3n":
				trip += [sep[0]]
			else:
				dipUnk += [sep[0]]
			line = infile.readline()
	
	# genotype triploids
	if len(trip) > 0:
		print("Genotyping triploids")
		subprocess.run(["mkdir", "tripGenos"])
		tripfiles = [x + ".bam" for x in trip] + [x + ".fastq" for x in trip]
		for f in tripfiles:
			subprocess.run(["mv", f, "tripGenos/"])
		# genotype triploids
		shellCommand = ["microhapWrap.py", "-t", str(threads), "-p", panelInfo[0], "-r", panelInfo[1], "-ploidy", "3", "-bam"]
		if len(panelInfo[3]) > 0:
			shellCommand += ["-pa", panelInfo[3]]
		subprocess.run(shellCommand, cwd = "./tripGenos/")
		
		# read in summary stats
		sumStatDict = {}
		with open("./tripGenos/summary_stats.txt", "r") as infile:
			line = infile.readline() # skip header
			line = infile.readline() # "Indiv" "GenotypeSuccess" "PropAlign" "ContamScore" "CountAlign" "TotalReads"
			while line:
				sep = line.rstrip().split("\t")
				for i in [1,3,4,5]:
					if sep[i] != "NA":
						sep[i] = float(sep[i])
				sumStatDict[sep[0]] = [sep[1], sep[3], sep[4], sep[5]] # gs, contam, align, total
				line = infile.readline()
		
		# make contam graph for triploids
		cScore = {}
		gS = {}
		allPeds = []
		for ind in sumStatDict:
			if re.match("rr?[0-9]NTC|initialNTC|f[0-9]NTC|qcNTC", ind):
				continue
			temp = sumStatDict[ind]
			if temp[1] == "NA" or temp[0] == "NA": # geno success should never be "NA", but just in case
				continue
			tempPed = re.sub("^initial|^f[0-9]|^rr?[0-9]|^qc|_[0-9]{4}$", "", ind)
			if tempPed not in allPeds:
				allPeds += [tempPed]
			cScore[tempPed] = cScore.get(tempPed, []) + [temp[1]]
			gS[tempPed] = gS.get(tempPed, []) + [temp[0]]
		
		colorPalette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
		fig, ax = plt.subplots()
		for i in range(0, len(allPeds)):
			ax.scatter(x = gS[allPeds[i]], y = cScore[allPeds[i]], c = colorPalette[i % len(colorPalette)], 
				label = allPeds[i], alpha = 0.4)
		lg = ax.legend(bbox_to_anchor = (1.05, -0.1), fancybox = True, framealpha = 0.4) # legend outside of plot
		ax.set_xlabel("Genotype success")
		ax.set_ylabel("Contamination score")
		plt.savefig("./tripGenos/triploid_ContamGraph.pdf", format = "pdf", bbox_extra_artists = (lg,), bbox_inches = "tight")

	# genotype (with full pipeline) diploids/unknown
	if len(glob.glob("*.bam")) > 0: # NTC may be left over but not in dipUnk
		print("Genotyping diploids and unknown ploidy samples")
		shellCommand = ["GTsEagle.py", "-panel", panel, "-t", str(threads), 
			"-s", str(minSuccess), "-pre", prefix, "--noAlign"]
		if bc is not None:
			shellCommand += ["-bc", bc]
		subprocess.run(shellCommand)

Main()
