#!/usr/bin/env python3
# wrapper to call genotypes with mtype2 (microTyper) and 
# calculate other related metrics important for genotyping
# also genotypes presence absence markers (e.g. sex markers)

import glob
import re
import sys
import os
import multiprocessing as mp
import subprocess

# uses wc to return number of reads in a 
#  fastq file (number of lines / 4)
def countReads(fastq):
	readCount = 0
	if os.path.exists(fastq):
		shellOut = subprocess.run(["wc", "-l", fastq], 
			stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True)
		readCount = int(shellOut.stdout.split(" ")[0]) / 4.0
	else:
		print("Warning: Did not find", fastq)
	return [fastq, readCount]

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
	
# uses samtools and cut to get alignment stats and
#  then count number of reads aligning to each locus
def countAlign(bam):
	alignCounts = {} # key is locus name, value is counts
	sumAlign = 0
	if os.path.exists(bam):
		shellOut = subprocess.run("samtools view " + bam + " | cut -f 3", 
			shell = True, stdout = subprocess.PIPE, universal_newlines = True)
		lines = shellOut.stdout.splitlines()
		for i in lines:
			if i == "*": # skip unmapped reads
				continue
			sumAlign += 1
			alignCounts[i] = alignCounts.get(i, 0) + 1
	return [bam, alignCounts, sumAlign]

# genotyping presence absence markers	
def callPAgenos(sampleName, totalAlign, totalAlignDict, paDict):
	# subtract any presence absence marker alignments from total alignment number
	for m in paDict:
		totalAlign -= totalAlignDict.get(m,0)
	# genotype
	output = []
	for m in paDict:
		paReads = totalAlignDict.get(m,0)
		tempDict = paDict[m]
		ratio = None
		geno = None
		if paReads == 0:
			ratio = "Inf"
			geno = tempDict[1]
		else:
			ratio = totalAlign / paReads
			if ratio <= tempDict[3]:
				geno = tempDict[0]
			elif ratio >= tempDict[4]:
				geno = tempDict[1]
			else:
				geno = "" # missing
		if totalAlign < tempDict[2]: # below min read count
			geno = ""
		output += [[sampleName, m, geno, str(ratio), str(paReads), str(totalAlign)]] # sample, locus, geno, ratio, paReads, alignedReads
	return output

# count reads in fastq file (single line sequences) that begin with a forward primer
def countFwd(fastq, fwdDict):
	# initialize with zeros so all are represented even if count is zeros
	fwdCountDict = {}
	for l in fwdDict:
		fwdCountDict[l] = 0
	with open(fastq, "r") as infile:
		line = infile.readline() # header
		while line:
			line = infile.readline() # read
			for l in fwdDict:
				if fwdDict[l].match(line):
					fwdCountDict[l] += 1
					break
			line = infile.readline() # spacer
			line = infile.readline() # quality
			line = infile.readline() # header
	return [re.sub("\.fastq$", "", fastq), fwdCountDict] # sampleName dictionary

def Main():
	
	# defaults
	fastqFiles = [] # -f
	bamFiles = [] # -bam
	align = True # whether to align the fastq files or not. Set to false if -b specified
	threads = None # -t, number of threads to use throughout, default is max available
	# inputs for mtype2
	pos = None # -p
	ref = None # -r 
	batch = 200 # -b 
	c = 0.99 # -c 
	eps = 0.01 # -eps 
	d = None # -d
	ploidy = 2 # -ploidy 
	paMarkInput = None # -pa, input file for pres/abs markers
	alignRef = None # -bt2ref, reference to pass to aligner
	outAlignByLocus = False # --alignByLocus, whether to output alignment counts by sample and locus
	fwd = None # -fwd, file with header line and "Locus \t fwdPrimerRegEx"

	# get command line inputs
	flag = 1
	while flag < len(sys.argv):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == "-f":
			flag += 1
			while flag < len(sys.argv):
				if sys.argv[flag][0] == "-":
					flag -= 1
					break
				fastqFiles += [sys.argv[flag]]
				flag += 1
		elif sys.argv[flag] == "-bam":
			align = False
			flag += 1
			while flag < len(sys.argv):
				if sys.argv[flag][0] == "-":
					flag -= 1
					break
				bamFiles += [sys.argv[flag]]
				flag += 1
		elif sys.argv[flag] == "-t":
			flag += 1
			threads = int(sys.argv[flag])
		elif sys.argv[flag] == "-p":
			flag += 1
			pos = sys.argv[flag]
		elif sys.argv[flag] == "-r":
			flag += 1
			ref = sys.argv[flag]
		elif sys.argv[flag] == "-b":
			flag += 1
			batch = int(sys.argv[flag])
		elif sys.argv[flag] == "-c":
			flag += 1
			c = float(sys.argv[flag])
		elif sys.argv[flag] == "-eps":
			flag += 1
			eps = float(sys.argv[flag])
		elif sys.argv[flag] == "-ploidy":
			flag += 1
			ploidy = int(sys.argv[flag])
		elif sys.argv[flag] == "-d":
			flag += 1
			d = int(sys.argv[flag])
		elif sys.argv[flag] == "-pa":
			flag += 1
			paMarkInput = sys.argv[flag]
		elif sys.argv[flag] == "-bt2ref":
			flag += 1
			alignRef = sys.argv[flag]
		elif sys.argv[flag] == "-fwd":
			flag += 1
			fwd = sys.argv[flag]
		elif sys.argv[flag] == "--alignByLocus":
			outAlignByLocus = True
		else:
			print("Error: option", sys.argv[flag], "not recognized.")
			return
		flag += 1
	# end while loop for command line arguments
	
	# default thread number
	if threads is None:
		try:
			threads = len(os.sched_getaffinity(0)) # -t, number of threads to use throughout, default is max available
		except:
			print("Error: could not detect number of threads available with your OS. Please specify number of threads using the -t argument")
			return
	
	if d is None:
		d = 5 * ploidy
	
	if pos is None:
		print("-p must be specified")
		return
	if ref is None:
		print("-r must be specified")
		return

	# get sample list if not input
	if align:
		if len(fastqFiles) == 0:
			fastqFiles = glob.glob("*.fastq")
			bamFiles = [re.sub("\.fastq$", ".bam", x) for x in fastqFiles]
	else:
		if len(bamFiles) == 0:
			bamFiles = glob.glob("*.bam")
		fastqFiles = [re.sub("\.bam$", ".fastq", x) for x in bamFiles]
	
	# create process pool
	print("Using ", str(threads), " thread(s).")
	process_pool = mp.Pool(processes=threads)
	
	# get total reads for each sample
	print("Counting total reads for each sample")
	totalReads = [process_pool.apply_async(countReads, args = [x]) for x in fastqFiles]
	
	totalReadsDict = {} # key is sample name value is total read count
	for r in totalReads:
		temp = r.get()
		totalReadsDict[re.sub("\.fastq$", "", temp[0])] = temp[1]
	del totalReads
	
	# align, if needed
	if align:
		print("Aligning reads")
		if alignRef is None:
			print("Error: Alignment reference prefix must be specified (-bt2ref).")
			return
		alignOut = [process_pool.apply_async(alignReads, args = [x, alignRef]) for x in fastqFiles]
		
		# output stderr (alignment stats)
		with open("alignment_output.txt", "w") as outfile:
			for i in alignOut:
				temp = i.get()
				outfile.write(temp[0] + "\n" + temp[1])

	# call genotypes (mtype2)
	print("Calling genotypes")
	subprocess.run(["mtype2", "-f"] + bamFiles + ["-p", pos, "-r", ref,
		"-o", "microhap_genotypes.txt", "-t", str(threads), "-b", str(batch), "-c", str(c), "-d", str(d),
		"-eps", str(eps), "-ploidy", str(ploidy)])
	
	# get alignment counts per locus
	print("Counting aligned reads per sample and locus")
	countsAligned = [process_pool.apply_async(countAlign, args = [x]) for x in bamFiles]
	alignReadsDict = {} # key is sample name value is dictionary with key locus name value count of aligned reads
	totalAlignDict = {} # key is sample name value is total number of aligned reads
	for r in countsAligned:
		temp = r.get()
		temp[0] = re.sub("\.bam$", "", temp[0])
		alignReadsDict[temp[0]] = temp[1]
		totalAlignDict[temp[0]] = temp[2]
	del countsAligned
	
	# genotype pres/abs markers
	paGenoSuccess = {} # keeping track for calculation of overall genotyping success
	if paMarkInput is not None:
		print("Genotyping presence/absence markers")
		paDict = {}
		numPA = 0
		# read in inputs
		with open(paMarkInput, "r") as infile:
			line = infile.readline() # skip header
			line = infile.readline() # Locus	PresGeno	AbsGeno	minReads	maxRatioPres	MinRatioAbs
			while line:
				sep = line.rstrip().split("\t")
				for i in range(3, 6, 1):
					sep[i] = float(sep[i])
				paDict[sep[0]] = sep[1:]
				numPA += 1
				line = infile.readline()
		# call genotypes
		paGenos = [process_pool.apply_async(callPAgenos, 
			args = [x, totalAlignDict[x], alignReadsDict[x], paDict]) for x in totalAlignDict]
		# write out
		with open("pres_abs_genotypes.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "Genotype", "Ratio", "LocusReads", "OtherAlignedReads"]) + "\n")
			for pa in paGenos:
				result = pa.get()
				genoSuccess = [0,numPA]
				sampleName = result[0][0]
				for line in result:
					if line[2] != "":
						genoSuccess[0] += 1
					outfile.write("\t".join(line) + "\n")
				paGenoSuccess[sampleName] = genoSuccess
		del paDict
		del paGenos
	
	# count reads beginnign with each fwd primer in each indiv
	if fwd is not None:
		fwdDict = {}
		print("Counting reads beginning with each forward primer")
		# read in primers
		with open(fwd, "r") as infile:
			line = infile.readline() # skip header
			line = infile.readline() # Locus	fwd
			while line:
				sep = line.rstrip().split("\t")
				fwdDict[sep[0]] = re.compile(sep[1])
				line = infile.readline()
		# count fwd primer reads
		fwdResults = [process_pool.apply_async(countFwd, args = [x, fwdDict]) for x in fastqFiles]
		# write out
		with open("fwd_primer_counts.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "FwdPrimerCounts"]) + "\n")
			for ind in fwdResults:
				result = ind.get()
				for loc in fwdDict:
					outfile.write("\t".join([result[0], loc, str(result[1][loc])]) + "\n")
		del fwdDict
		del fwdResults
		
	
	process_pool.close()  # close the pool to new tasks
	process_pool.join()  # join the processes
	
	# calc summary stats
	# genotyping success, proportion reads align, contamination detection
	print("Calculating summary statistics")
	genoSuccess = {}
	contam = {}
	for s in totalAlignDict: # inititating dictionaries
		genoSuccess[s] = [0,0] # [num successful, total number of loci]
		contam[s] = 0 # sum contamination metric, note that number of loci is num successful in line above
	with open("microhap_genotypes.txt", "r") as infile:
		line = infile.readline() # skip header
		line = infile.readline() # Sample.bam	Locus	Allele1	...	Allele1_count	...	p
		while line:
			sep = line.split("\t")
			sn = re.sub("\.bam$", "", sep[0]) # sample name
			genoSuccess[sn][1] += 1
			if sep[2] != "":
				genoSuccess[sn][0] += 1
				# contamination score
				#  calculated as reads of |(Allele1 / total aligned reads) - (Allele1 dosage / ploidy)|
				a1 = sep[2]
				a1Dose = 1.0
				for i in range(1, ploidy, 1):
					if sep[2 + i] == a1:
						a1Dose += 1
				contam[sn] += abs((float(sep[ploidy + 2]) / alignReadsDict[sn][sep[1]]) - (a1Dose / ploidy))
			line = infile.readline()
	
	with open("summary_stats.txt", "w") as outfile:
		outfile.write("\t".join(["Indiv", "GenotypeSuccess", "PropAlign", "ContamScore", "CountAlign", "TotalReads"]) + "\n")
		for s in totalAlignDict: # for each sample
			gs = genoSuccess[s]
			paGS = paGenoSuccess.get(s, [0,0])
			propAlign = "NA"
			if totalReadsDict[s] > 0:
				propAlign = totalAlignDict[s] / totalReadsDict[s]
			co = "NA"
			if gs[0] > 0:
				co = contam[s] / gs[0]
			outfile.write("\t".join([s, str((gs[0] + paGS[0])/(gs[1] + paGS[1])), str(propAlign), str(co), str(totalAlignDict[s]), str(totalReadsDict[s])]) + "\n")
	
	# optionally output count of aligned reads per sample and locus
	if outAlignByLocus:
		with open("align_by_ind_locus.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "CountAlign"]) + "\n")
			for s in alignReadsDict:
				temp = alignReadsDict[s]
				for l in temp:
					outfile.write("\t".join([s, l, str(temp[l])]) + "\n")

Main()
