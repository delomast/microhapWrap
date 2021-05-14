#!/usr/bin/env python3
# wrapper for processing panels with microTyper and microhapWrap
# for the Eagle Fish Genetics Laboratory

import re
import sys
import os
import glob
import subprocess
import matplotlib.pyplot as plt

def Main():
	# defaults
	panel = None # -panel, panel name
	threads = len(os.sched_getaffinity(0)) # -t, number of threads to use throughout, default is max available
	bc = None # -bc, bcSplitFile.csv
	minSuccess = 0.9 # -s, minimum proportion of genotypes called to be "successful"
	prefix = "GTsEagle" # -pre, prefix for output files
	
	# wells as numbered in the bcSplit file from 1 - 96
	wellPosE = ["A01", "B01", "C01", "D01", "E01", "F01", "G01", "H01", "A02", "B02", "C02", "D02", "E02", 
		"F02", "G02", "H02", "A03", "B03", "C03", "D03", "E03", "F03", "G03", "H03", "A04", "B04", "C04", "D04", 
		"E04", "F04", "G04", "H04", "A05", "B05", "C05", "D05", "E05", "F05", "G05", "H05", "A06", "B06", "C06", 
		"D06", "E06", "F06", "G06", "H06", "A07", "B07", "C07", "D07", "E07", "F07", "G07", "H07", "A08", "B08", 
		"C08", "D08", "E08", "F08", "G08", "H08", "A09", "B09", "C09", "D09", "E09", "F09", "G09", "H09", "A10", 
		"B10", "C10", "D10", "E10", "F10", "G10", "H10", "A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11", 
		"A12", "B12", "C12", "D12", "E12", "F12", "G12", "H12"]

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
			"/home/efglserv/software/GTseq/microhapInputs/presAbsInputs/Sna302_pres_abs.txt"],
			
		"Svi151" : ["/home/efglserv/software/GTseq/microhapInputs/posFiles/Svi151_pos_final.txt", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/Svi151/wall_1.fasta", 
			"/home/efglserv/software/GTseq/microhapInputs/refSeqs/Svi151/wall_1", 
			""] # no pres abs markers for Svi151
	}
		
	if panel not in panelDict:
		print("Valid panel names are:")
		for p in panelDict:
			print(p)
		return
		
	panelInfo = panelDict[panel]
	
	# call microhapWrap
	shellCommand = ["microhapWrap.py", "-t", str(threads), "-p", panelInfo[0], "-r", panelInfo[1], 
		"-bt2ref", panelInfo[2]]
	if len(panelInfo[3]) > 0:
		shellCommand += ["-pa", panelInfo[3]]
	subprocess.run(shellCommand)
	
	# now organize and create lab specific outputs
	# different genotype files
	iGenos = {} # key as individual name and value as [[locus1, a1, a2], [locus2, a1, a2], ...]
	fGenos = {}
	qcGenos = {}
	rrGenos = {}
	otherGenos = {}
	
	locusSuccess = {} # locus name : successful
	locusTotal = {} # locus name : total

	with open("microhap_genotypes.txt", "r") as infile:
		line = infile.readline() # skip header
		line = infile.readline() # Sample.bam	Locus	Allele1	...	Allele1_count	...	p
		while line:
			sep = line.split("\t")
			sep[0] = re.sub("\.bam$", "", sep[0]) # sample name
			if sep[2] == "": 
				sep[2] = "000"
				sep[3] = "000"
			else:
				locusSuccess[sep[1]] = locusSuccess.get(sep[1], 0) + 1
			if re.match("initial", sep[0]):
				locusTotal[sep[1]] = locusTotal.get(sep[1], 0) + 1
				if sep[2] != "000":
					locusSuccess[sep[1]] = locusSuccess.get(sep[1], 0) + 1
				iGenos[sep[0]] = iGenos.get(sep[0], []) + [sep[1:4]]
			elif re.match("f[0-9]", sep[0]):
				fGenos[sep[0]] = fGenos.get(sep[0], []) + [sep[1:4]]
			elif re.match("qc", sep[0]):
				qcGenos[sep[0]] = qcGenos.get(sep[0], []) + [sep[1:4]]
			elif re.match("rr?[0-9]", sep[0]):
				rrGenos[sep[0]] = rrGenos.get(sep[0], []) + [sep[1:4]]
			else:
				otherGenos[sep[0]] = otherGenos.get(sep[0], []) + [sep[1:4]]
			line = infile.readline()
	
	# read in summary stats
	sumStatDict = {}
	with open("summary_stats.txt", "r") as infile:
		line = infile.readline() # skip header
		line = infile.readline() # "Indiv" "GenotypeSuccess" "PropAlign" "ContamScore" "CountAlign" "TotalReads"
		while line:
			sep = line.rstrip().split("\t")
			for i in [1,3,4,5]:
				if sep[i] != "NA":
					sep[i] = float(sep[i])
			sumStatDict[sep[0]] = [sep[1], sep[3], sep[4], sep[5]] # gs, contam, align, total
			line = infile.readline()
	
	# genotype files
	if len(iGenos) > 0:
		with open(prefix + "_initialMicrohaps.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "Allele1", "Allele2"]) + "\n")
			for ind in iGenos:
				for line in iGenos[ind]:
					outfile.write("\t".join([re.sub("^initial", "", ind)] + line) + "\n")
	if len(fGenos) > 0:
		with open(prefix + "_ALLfillin_Microhaps.txt", "w") as outfile1:
			with open(prefix + "_successFillin_Microhaps.txt", "w") as outfile2:
				outfile1.write("\t".join(["Indiv", "Locus", "Allele1", "Allele2"]) + "\n")
				outfile2.write("\t".join(["Indiv", "Locus", "Allele1", "Allele2"]) + "\n")
				for ind in fGenos:
					for line in fGenos[ind]:
						outfile1.write("\t".join([re.sub("^f[0-9]", "", ind)] + line) + "\n")
					if sumStatDict[ind][0] >= minSuccess:
						for line in fGenos[ind]:
							outfile2.write("\t".join([re.sub("^f[0-9]", "", ind)] + line) + "\n")
	if len(qcGenos) > 0:
		with open(prefix + "_qcMicrohaps.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "Allele1", "Allele2"]) + "\n")
			for ind in qcGenos:
				for line in qcGenos[ind]:
					outfile.write("\t".join([re.sub("^qc", "", ind)] + line) + "\n")
	if len(rrGenos) > 0:
		with open(prefix + "_rrMicrohaps.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "Allele1", "Allele2"]) + "\n")
			for ind in rrGenos:
				for line in rrGenos[ind]:
					outfile.write("\t".join([re.sub("^rr?[0-9]", "", ind)] + line) + "\n")
	if len(otherGenos) > 0:
		with open(prefix + "_otherMicrohaps.txt", "w") as outfile:
			outfile.write("\t".join(["Indiv", "Locus", "Allele1", "Allele2"]) + "\n")
			for ind in otherGenos:
				for line in otherGenos[ind]:
					outfile.write("\t".join([ind] + line) + "\n")
	
	# poor markers
	if len(locusTotal) > 0:
		with open(prefix + "_poorMarkers.csv", "w") as outfile:
			outfile.write("Markers with under 80% success in the initial samples,NoCall Rate\n")
			for locus in locusTotal:
				prop = 1 - locusSuccess.get(locus, 0) / locusTotal[locus]
				if prop > 0.2:
					outfile.write(",".join([locus, str(prop)]) + "\n")
	
	# make contam graph
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
	ax.legend(framealpha = 0)
	ax.set_xlabel("Genotype success")
	ax.set_ylabel("Contamination score")
	plt.savefig(prefix + "_ContamGraph.pdf", format = "pdf")
	
	
	# read in barcode file
	if bc is None:
		bc = glob.glob("BC*.csv") + glob.glob("bc*.csv") + glob.glob("barcode*.csv")
		if len(bc) < 1:
			print("Error: no bcSplit file specified or detected.")
			return
		bc = bc[0]
		print("Using auto-detected bcsplit file", bc)
	
	bcDict = {} # key is type + Name (e.g., initialOtsSAWT21S_0001), value is barcode fields in order
	allTrays = []
	with open(bc, "r") as infile:
		line = infile.readline() # skip header
		line = infile.readline() # SampleID	SampleType	SampleStatus	PlateID	i7_name	i7_sequence	i5_name	i5_sequence
		while line:
			sep = line.rstrip().split(",")
			if len(sep) < 8 or len(sep[0]) < 1: # skip blank or incorrectly formatted lines
				continue
			if sep[3] not in allTrays:
				allTrays += [sep[3]]
			sep[6] = int(sep[6])
			bcDict[sep[1] + sep[0]] = sep
			line = infile.readline()

	# dGaps
	with open(prefix + "_dGaps.csv", "w") as outfile:
		outfile.write(",".join(["Indiv", "sampleType", "sampleStatus", "propSuccess", "PlateID", "Well"]) + "\n")
		for ind in sumStatDict:
			if re.match("rr?[0-9]NTC|initialNTC|f[0-9]NTC|qcNTC", ind):
				continue
			if sumStatDict[ind][0] < minSuccess:
				temp = bcDict[ind]
				outfile.write(",".join([temp[0], temp[1], temp[2], str(sumStatDict[ind][0]), temp[3], wellPosE[temp[6] - 1]]) + "\n")

	# dGapsMap
	allTrayFails = {}
	letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
	with open(prefix + "_dGapsMap.csv", "w") as outfile:
		for tray in allTrays:
			tempFails = {}
			# this is not very efficient, but it works and isn't too slow given sample numbers
			for ind in sumStatDict:
				if re.match("rr?[0-9]NTC|initialNTC|f[0-9]NTC|qcNTC", ind):
					continue
				if bcDict[ind][3] == tray and sumStatDict[ind][0] < minSuccess:
					tempFails[bcDict[ind][6]] = ind
					allTrayFails[bcDict[ind][6]] = allTrayFails.get(bcDict[ind][6], 0) + 1
			outfile.write(tray + ",1,2,3,4,5,6,7,8,9,10,11,12\n")
			for i in range(0, 8, 1):
				line = [letters[i]]
				for j in range((12 * i) + 1, (12 * (i + 1)) + 1, 1):
					line += [tempFails.get(j, "")]
				outfile.write(",".join(line) + "\n")
			outfile.write("\n")
		outfile.write("AllTrays,1,2,3,4,5,6,7,8,9,10,11,12\n")
		for i in range(0, 8, 1):
			line = [letters[i]]
			for j in range((12 * i) + 1, (12 * (i + 1)) + 1, 1):
				line += [str(allTrayFails.get(j, ""))]
			outfile.write(",".join(line) + "\n")
	
	# stats by tray (report)
	with open(prefix + "_traySummary.csv", "w") as outfile:
		outfile.write("Tray summary excluding NTCs\n")
		outfile.write("PlateID,totalReads,alignedReads,meanPropAlign,n,numFailed,propSamplesSuccess\n")
		for tray in allTrays:
			totalReads = 0
			alignedReads = 0
			propAlign = 0
			nPropalign = 0
			n = 0
			numFail = 0
			# this is not very efficient, but it works and isn't too slow given sample numbers
			for ind in sumStatDict: # gs, contam, align, total
				if re.match("rr?[0-9]NTC|initialNTC|f[0-9]NTC|qcNTC", ind):
					continue
				if bcDict[ind][3] == tray:
					temp = sumStatDict[ind]
					totalReads += temp[3]
					alignedReads += temp[2]
					if temp[3] > 0:
						propAlign += temp[2] / temp[3]
						nPropalign += 1
					n += 1
					if sumStatDict[ind][0] < minSuccess:
						numFail += 1
			if nPropalign > 0:
				propAlign = propAlign / nPropalign
			else:
				propAlign = "NA"
			if n > 0:
				propSampSuccess = (n - numFail) / nPropalign
			else:
				propSampSuccess = "NA"
			outfile.write(",".join([tray, str(totalReads), str(alignedReads), str(propAlign), str(n), str(numFail), str(propSampSuccess)]) + "\n")
		outfile.write("\nNTC summary\nIndiv,propSuccess\n")
		for ind in sumStatDict: # gs, contam, align, total
			if re.match("rr?[0-9]NTC|initialNTC|f[0-9]NTC|qcNTC", ind):
				outfile.write(",".join([ind, str(sumStatDict[ind][0])]) + "\n")
		
	# Dfunc
	# separate from other geno output files b/c want geno output files to be produced even without
	# barcode file present
	with open(prefix + "_DfuncItUp.txt", "w") as outfile:
		outfile.write("\t".join(["Indiv", "Type", "Status", "Plate", "Well", "Locus", "Allele1", "Allele2"]) + "\n")
		allGenos = {**iGenos, **fGenos, **qcGenos, **rrGenos, **otherGenos} # combining all dictionaries
		for ind in allGenos:
			for line in allGenos[ind]:
				outfile.write("\t".join(bcDict[ind][0:4] + [wellPosE[bcDict[ind][6] - 1]] + line) + "\n")
	
	print("Kaw! Kaw!")

Main()
