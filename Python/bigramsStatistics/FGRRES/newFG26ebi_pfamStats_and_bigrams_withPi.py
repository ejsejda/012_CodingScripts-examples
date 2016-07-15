#!/usr/bin/python
#PIPE SEPARATED PAIRS IN BIGRAMS
##################################################################
# USED PYTHON LIBRARY
##################################################################
from collections import OrderedDict
import sys


#################################################################
#  INPUT AND OUTPUT FILES - PLEASE ADJUST THE INPUT DIRECTORY AND THE OUTPUT TO YOUR NEEDS
#################################################################
workDir1 = "/home/ela/Project/pfam(sinceJune2013)" 

infile1 = "%s/input/solveOverlapping/goldFiles/newFusariumGraminearumEBI26_OverlappingSolved_GOLD.txt" % (workDir1)#raw file of pfam domain in newFG26ebi in order

infile2 = "%s/input/pfam27.0/pfamA.txt" % (workDir1)# raw file taken from PFAM 27.0 database

outfile1 = "%s/output/bigramsNewFG26ebi/pfam_to_Duf_map_forAllPfam27.0_pi.txt" % (workDir1)
outfile2 = "%s/output/bigramsNewFG26ebi/pfamDomCount_inNewFG26ebi_withUniquePfamNo_pi.txt" % (workDir1)
outfile3 = "%s/output/bigramsNewFG26ebi/pfamDomCount_inNewFG26ebi_with_Only_1_UniquePfamNo_pi.txt" % (workDir1)
outfile4 = "%s/output/bigramsNewFG26ebi/pfamDomCount_inNewFG26ebi_withOnlyOnePfamDom_pi.txt" % (workDir1)
outfile5 = "%s/output/bigramsNewFG26ebi/pfamDomStatisticNewFG26ebi_pi.txt" % (workDir1)
outfile6 = "%s/output/bigramsNewFG26ebi/uniquePfamInNewFG26ebi_withPfamDescr_withDufs_pi.txt" % (workDir1)
outfile7 = "%s/output/bigramsNewFG26ebi/pfamDomBigrams_domsOrderNoMatter/all_UniqueBigramsOrderNoMatter_pi.txt" % (workDir1)
outfile8 = "%s/output/bigramsNewFG26ebi/pfamDomBigrams_domsOrderNoMatter/homo_UniqueBigramsOrderNoMatter_pi.txt" % (workDir1)
outfile9 = "%s/output/bigramsNewFG26ebi/pfamDomBigrams_domsOrderNoMatter/hetero_UniqueBigramsOrderNoMatter_pi.txt" % (workDir1)
outfile10 = "%s/output/bigramsNewFG26ebi/pfamDomBigrams_domsOrderMatters/all_UniqueBigrams_OrderMatters_pi.txt" % (workDir1)
outfile11 = "%s/output/bigramsNewFG26ebi/pfamDomBigrams_domsOrderMatters/homo_UniqueBigrams_OrderMatters_pi.txt" % (workDir1)
outfile12 = "%s/output/bigramsNewFG26ebi/pfamDomBigrams_domsOrderMatters/hetero_UniqueBigrams_OrderMatters_pi.txt" % (workDir1)

######################################################################################
#          OPENNING FILES FOR READING
####################################################################################
fh1 = open (infile1, "r") 
fh2 = open (infile2, "r")
fh_out1 = open (outfile1, "w")
fh_out2 = open (outfile2, "w")
fh3 = open (outfile2, "r")
fh_out3 = open (outfile3, "w")
fh_out4 = open (outfile4, "w")
fh_out5 = open (outfile5, "w")
fh_out6 = open (outfile6, "w")
fh_out7 = open (outfile7, "w")
fh_out8 = open (outfile8, "w")
fh_out9 = open (outfile9, "w")
fh_out10 = open (outfile10, "w")
fh_out11 = open (outfile11, "w")
fh_out12 = open (outfile12, "w")

######################################################################################
# DECLARATION OF GLOBAL DICTIONARIES AND SETS
######################################################################################
pfam_to_Duf = dict()

duf_count = set()#set of all DUFs in newFG26ebi proteome and the whole set includes all DUFs

pfamId_Desc = dict() #dictionary where keys are PFAM or DUF (where present) domains and values are the description of the domain

newFG26ebi_prot = dict()#dictionary where keys are newFG26ebi_Id and value is a list of total pfam domain within given newFG26ebi protein, there are ....newFG26ebi proteins that have PFAM(including DUF  where exist) domains

newFG26ebi_uniqPfam = dict()#dictionary where keys are newFG26ebi_Id and value is a list of unique pfam domain within given newFG26ebi protein

pfamIds = dict()#dictionary where keys are pfam/or DUF ids and the value are total number of newFG26ebi_Ids where the given PFAM/or DUF is 	present, so value can be several the same newFG26ebi_ids

#pfamStartStop = dict()
OneDufDomnewFG26ebi= set()
uniqDufs_InOneDufDomnewFG26ebi = set()
pfamNum = dict()

################################################################################################
#1. MAPPING PFAM 27.0 ACCESION NUMBERS TO DUF ID FROM PFAM27.0 DATABASE
#	Creating dictionary for PfamIds in pfam 27.0 to Duf_ID
#2. CREATING FILE WITH MAPPING DICTIONARY OF PFAM IDs AND DUFs IDs FOR WHOLE PFAM VERSION 27.0
#################################################################################################
def creatingPfamDufsMap_file():
	
	for line in fh2.readlines():
		line = line.split('\t')
		pfam_Acc = line[1].strip()
		pfam_Desc = line[2].strip()
		if line[2].startswith('DUF'):
			pfam_to_Duf.setdefault(pfam_Acc,pfam_Desc)

	for key, value in pfam_to_Duf.items():
	#fh_out1.write("%s\t%s\n" %(key,str(value).replace('[\'','').replace('\']','').replace('\'','').replace(',','\t').strip()))
		fh_out1.write("%s\t%s\n" %(key,value))
	fh_out1.close()
	return (pfam_to_Duf)
##################################################################################################	
#READING THE ROW FILE WITH newFG26ebi_Ids WITH PFAM Acc FROM PFAM 27.
#THIS IS A RAW INPUT FILE GENERATED AFTER SOLVING THE OVERLAPPING OF DOMAINS IN THE newFG26ebi PROTEINS
# generating following dictionaries:
#1. pfamId_Desc, 2. newFG26ebi_prot, 3.newFG26ebi_uniqPfam, 4.pfamIds and set: duf_count
##################################################################################################
def readingGoldFile():

	for line in fh1.readlines():
		line = line.split('\t')
		newFG26ebiId = line[0].strip()#first column in the file1 is newFG26ebiSG_id
		pfamAcc = line[1][:7].strip()#second column in the file1 is pfam accesion number and here we remove everything after 7 						      character of pfam accesion no
		if pfamAcc in pfam_to_Duf:
			pfamAcc = pfam_to_Duf[pfamAcc]#converting pfamId to DUF id if present
			duf_count.add(pfamAcc)
		pfamDesc = line[2].strip()
		pfamStart = line[3][:-2].strip()
		pfamEnd = line[4][:-3].strip()
		pfamInfo = str(pfamAcc) +"_" + str(pfamStart) + "-" + str(pfamEnd)#concatenating start and end domain into pfamAccesion

		pfamId_Desc.setdefault(pfamAcc, pfamDesc)
		newFG26ebi_prot.setdefault(newFG26ebiId,list())
		newFG26ebi_uniqPfam.setdefault(newFG26ebiId,set())
		#pfamStartStop.setdefault(newFG26ebiId, list())#setting newFG26ebiId as a key in pfam dictionary and value as a list as we take into account  every pfam domain per newFG26ebi protein. Thus, the same pfam accesion is allowed in one newFG26ebi protein.
		newFG26ebi_prot[newFG26ebiId].append(pfamAcc)
		newFG26ebi_uniqPfam[newFG26ebiId].add(pfamAcc)
		#pfamStartStop[newFG26ebiId].append(pfamInfo)#filling the value - list with pfam accesion number for each newFG26ebi protein.
		pfamIds.setdefault(pfamAcc,list())
		pfamIds[pfamAcc].append(newFG26ebiId)
	return(pfamId_Desc, newFG26ebi_prot, newFG26ebi_uniqPfam, pfamIds, duf_count)
##########################################################################################################################################
#####################################################################################################
#GENERATING FILE: newFG26ebiSG proteins with PFAM domains in order per sequence including total 
#and unique numbers of PFAM (DUF) domains per newFG26ebi protein. It returns following sets: 
#1.OneDufDomnewFG26ebi, 2.uniqDufs_InOneDufDomnewFG26ebi, and dictionary with statistic: pfamNum
#####################################################################################################
def pfamDomainsInnewFG26ebi():
	print "******************************************************************************"


	newFG26ebi_prot_sorted_by_value =OrderedDict(sorted(newFG26ebi_prot.items(), key = lambda x: len(x[1]), reverse= True))

	for key, value in newFG26ebi_prot_sorted_by_value.items():
		if key in newFG26ebi_uniqPfam.keys():
			#print key,'\t', len(value),'\t', str(value).replace('[\'','').replace('\']','').replace('\'','').replace(',','\t')
			fh_out2.write("%s\t%s\t%s\t%s\n" %(key,len(value),len(newFG26ebi_uniqPfam[key]),str(value).replace('[\'','').replace('\']','').replace('\'','').replace(', ','\t').strip()))
			if len(value) > 1 and len(newFG26ebi_uniqPfam[key])== 1:
				fh_out3.write("%s\t%s\t%s\t%s\n" %(key,len(value),len(newFG26ebi_uniqPfam[key]),str(value).replace('[\'','').replace('\']','').replace('\'','').replace(', ','\t').strip()))

			if len(value) == 1 and len(newFG26ebi_uniqPfam[key])== 1:
				fh_out4.write("%s\t%s\t%s\t%s\n" %(key,len(value),len(newFG26ebi_uniqPfam[key]),str(value).replace('[\'','').replace('\']','').replace('\'','').replace(', ','\t').strip()))
			
			if len(value) == 1 and str(value).replace('[\'','').replace('\']','').startswith("DUF"):
				if str(value) not in uniqDufs_InOneDufDomnewFG26ebi:
					uniqDufs_InOneDufDomnewFG26ebi.add(str(value))
				if key not in OneDufDomnewFG26ebi:
					OneDufDomnewFG26ebi.add(key)
			
			prot = key
			l = len(value)
			pfamNum.setdefault(l,set()) #l(number of pfam domains in given newFG26ebi protein) is a key in pfamNum dictionary, where value is a set(or can be list in this situation) of newFG26ebi ids(proteins) having those number of pfam within.
			pfamNum[l].add(prot)


	print "Unique Duf domains in all Single duf newFG26ebi proteins:", len(uniqDufs_InOneDufDomnewFG26ebi)
	print "No of one DUFs domain unique newFG26ebi proteins:", len(OneDufDomnewFG26ebi)
		
	for k, v in pfamNum.items():
		#print k, len(v)
		fh_out5.write("%s\t%s\n" %(k, len(v)))		

	fh_out2.close()
	fh_out3.close()
	fh_out4.close()		
	fh_out5.close()

	return(OneDufDomnewFG26ebi, uniqDufs_InOneDufDomnewFG26ebi, pfamNum)
############################################################################################################################################################################################################################################################################
# GENERATING FILE: List of unique PFAM(or DUFs) sorted by their number in newFG26ebi proteome,their function and list of newFG26ebiIds where they are present###
##################################################################################################################################################
def uniquePfamDufDomainsInnewFG26ebi():
	pfamIds_sorted_by_value =OrderedDict(sorted(pfamIds.items(), key = lambda x: len(x[1]), reverse= True))
	for key, value in pfamIds_sorted_by_value.items():
		#if key in pfamId_Desc.items():
		print key, len(value), pfamId_Desc[key]
		fh_out6.write("%s\t%s\t%s\t%s\n" %(key, len(value), pfamId_Desc[key], str(value).replace('[\'','').replace('\']','').replace('\'','').replace(', ','\t').strip()))
	print "Total number of unique PFAM domains within newFG26ebi proteome is:", len(pfamIds)
	fh_out6.close()
####################################################################################################

creatingPfamDufsMap_file()
readingGoldFile()
pfamDomainsInnewFG26ebi()
uniquePfamDufDomainsInnewFG26ebi()

####################################################################################################
#  CREATING ALL PFAM/(DUF) BIGRAMS PRESENT IN newFG26ebi PROTEOME###########################################
####################################################################################################
sortedBigrams = dict()
sortedBigrams_inUnqProt = dict()
allBigrams = dict()
allBigrams_inUnqProt = dict()
bigramCount = 0

hetero_pairs= dict()
hetero_pairs_inUniqProt = dict()
all_hetero_pairs = dict()
all_hetero_pairs_inUniqProt = dict()

homo_pairs = dict()
homo_pairs_inUniqProt = dict()
all_homo_pairs = dict()
all_homo_pairs_inUniqProt = dict()

hetero_bigramCount = 0
homo_bigramCount = 0

orderNoMatterBigrams=set()
orderMatterBigrams = set()

for l in fh3.readlines():
	line = l.strip().split('\t')

	for i in range(3, len(line)-1):
		a = line[i].strip()
		b = line[i+1].strip()
		bigramCount = bigramCount+1
		#bigram = a+"|"+b
		pair = [a,b]
	
		#str_pair = str(pair).replace('[\'','').replace('\', \'','\t').replace('\']','')
		pi_pair = str(pair).replace('[\'','').replace('\', \'','|').replace('\']','')
		#print pi_pair
		sortedPair = str(sorted(pi_pair)).replace('[\'','').replace('\', \'','|').replace('\']','')
		#print str_pair, sortedPair
		
		if pi_pair not in orderMatterBigrams:
			orderMatterBigrams.add(pi_pair)
			
		if sortedPair not in orderNoMatterBigrams:
			orderNoMatterBigrams.add(sortedPair)

		sortedBigrams.setdefault(sortedPair, list()).append(line[0])
		sortedBigrams_inUnqProt.setdefault(sortedPair, set()).add(line[0])

		allBigrams.setdefault(pi_pair, list()).append(line[0])
		allBigrams_inUnqProt.setdefault(pi_pair, set()).add(line[0])
		
		if a == b:
			homo_pairs.setdefault(sortedPair, list()).append(line[0])
			homo_pairs_inUniqProt.setdefault(sortedPair, set()).add(line[0])
			all_homo_pairs.setdefault(pi_pair, list()).append(line[0])
			all_homo_pairs_inUniqProt.setdefault(pi_pair, set()).add(line[0])
		else:
			hetero_pairs.setdefault(sortedPair, list()).append(line[0])
			hetero_pairs_inUniqProt.setdefault(sortedPair, set()).add(line[0])
			all_hetero_pairs.setdefault(pi_pair, list()).append(line[0])
			all_hetero_pairs_inUniqProt.setdefault(pi_pair, set()).add(line[0])

print "No of unique bigrams where order matters:", len(orderMatterBigrams)
print "No of unique bigrams where order doesn't matter:", len(orderNoMatterBigrams)
print len(sortedBigrams)
print "Number of all unique bigrams:",len(sortedBigrams)
print "Number of unique hetero bigrams:",len(hetero_pairs)
print "Number of unique homo bigrams:",len(homo_pairs) 

##################################################################################################################
def bigramsStatistic(bigrams1,bigrams2,outputFile):
	bigrams_sorted_by_value_withAllnewFG26ebi = OrderedDict(sorted(bigrams1.items(), key = lambda x: len(x[1]), reverse= True))
	bigrams_sorted_by_value_withUniqnewFG26ebi = OrderedDict(sorted(bigrams2.items(), key = lambda x: len(x[1]), reverse= True))
	for key, value in bigrams_sorted_by_value_withAllnewFG26ebi.items():
		if key in bigrams_sorted_by_value_withUniqnewFG26ebi.keys():
			#print key, len(value), len(bigrams_sorted_by_value_withUniqnewFG26ebi[key]), str(value).replace('[\'','').replace('\', \'','\t').replace('\']','')
			
			outputFile.write("%s\t%s\t%s\t%s\n" %(key, len(value), len(bigrams_sorted_by_value_withUniqnewFG26ebi[key]), str(bigrams_sorted_by_value_withUniqnewFG26ebi[key]).replace('set([\'','').replace('\', \'','\t').replace('\'])','')))
	outputFile.close()
####################################################################################################################
bigramsStatistic(sortedBigrams,sortedBigrams_inUnqProt,fh_out7)#all_bigrams where order of domains in bigrams doesn't matter
bigramsStatistic(homo_pairs,homo_pairs_inUniqProt,fh_out8)#homo_bigrams where order of domains in bigrams doesn't matter
bigramsStatistic(hetero_pairs,hetero_pairs_inUniqProt,fh_out9)#hetero_bigrams where order of domains in bigrams doesn't matter
bigramsStatistic(allBigrams,allBigrams_inUnqProt,fh_out10) #all_bigrams where order of domains in bigrams matters
bigramsStatistic(all_homo_pairs,all_homo_pairs_inUniqProt,fh_out11)#homo_bigrams where order of domains in bigrams matters
bigramsStatistic(all_hetero_pairs,all_hetero_pairs_inUniqProt,fh_out12)#hetero_bigrams where order of domains in bigrams matters
