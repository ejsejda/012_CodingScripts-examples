#!/usr/bin/python

from collections import OrderedDict
import sys
import collections

#Please adjust input directory files for your needs:
workDir1 = "/home/ela/Project/pfam(sinceJune2013)" 

infile_pfam27 = "%s/input/pfam27.0/pfamA.txt"% (workDir1)
infile1 = "%s/input/taxonomyInFG_nov2014/001-results_004_testResults_Fixed-finished-14Nov.txt"% (workDir1)
infile2 = "%s/input/taxonomyInFG_nov2014/all_no_pathogens.txt"% (workDir1)
infile3 = "%s/input/taxonomyInFG_nov2014/fungi_pathogen.txt"% (workDir1)
infile4 = "%s/input/taxonomyInFG_nov2014/plant_other_symbionts.txt"% (workDir1)
infile5 = "%s/input/taxonomyInFG_nov2014/plant_pathogens.txt"% (workDir1)
infile6 = "%s/input/taxonomyInFG_nov2014/animal_pathogens.txt"% (workDir1)



outfile1 = "%s/output/pfam27_taxa_Nov2014_copyJan2015/pfam27ToDUFs_copyJan2015.txt" % (workDir1)
outfile2 = "%s/output/pfam27_taxa_Nov2014_copyJan2015/DUF_Ids_within_unclassified_copyJan2015.txt" % (workDir1)
outfile3 = "%s/output/pfam27_taxa_Nov2014_copyJan2015/all_Pfam_Ids_within_unclassified_copyJan2015.txt" % (workDir1)
outfile_3 = "%s/output/pfam27_taxa_Nov2014_copyJan2015/dufIds_InFungi_and_PlantaOnly_copyJan2015.txt" % (workDir1)
outfile_4 = "%s/output/pfam27_taxa_Nov2014_copyJan2015/dufIds_InFungi_and_PlantaOnly_WithSpecies_copyJan2015.txt" % (workDir1)
outfile4 = "%s/output/pfam27_taxa_Nov2014_copyJan2015/dufId_in_fungiOnly_copyJan2015.txt" % (workDir1)
fh1 = open (infile1, "r")
fhOut2 = open(outfile2, "a")
fhOut3 = open(outfile3, "a")
fhOut_3 = open (outfile_3, "a")
fhOut_4 = open (outfile_4, "a")
fhOut4 = open (outfile4, "w")

################################################################################################
#1. MAPPING PFAM ACCESION NUMBERS TO DUF ID FROM PFAM DATABASE
#	Creating dictionary for PfamIds in pfam  to Duf_ID
#2. CREATING FILE WITH MAPPING DICTIONARY OF PFAM IDs AND DUFs IDs FOR WHOLE PFAM GIVEN VERSION
#################################################################################################
def creatingPfamDufsMap_file(infile, outfile, mappingDict):
	print " *************************************************"
	print "  Start mapping "
	fhin = open (infile, "r")
	fhout = open (outfile, "w")
	mappingDic = dict()
	for line in fhin.readlines():
		line = line.split('\t')
		pfam_Acc = line[1].replace('\'','').strip()
		pfam_Desc = line[2].replace('\'','').strip()
		if pfam_Desc.startswith('DUF'):
			mappingDic.setdefault(pfam_Acc,pfam_Desc)
			mappingDict[pfam_Acc] = pfam_Desc

	for key, value in mappingDic.items():
		fhout.write("%s\t%s\n" %(key,str(value).replace('[\'','').replace('\']','').replace('\'','').replace(',','\t').strip()))
		#fhout.write("%s\t%s\n" %(key,value))
		#print key, value
	print "there are %s pfam IDs in Pfam 27 database converted to DUF ids" %len(mappingDict)
	fhout.close()
	fhin.close()
	return (mappingDic)
#########################################################################################################
#########################################################################################################
#########################################################################################################
#declaring empty dictionaries
pfam27_to_Duf = dict()
dufFungiDic = dict()
allFungiDic = dict()
fungiWithDuf = dict()#keys are fungi name and values are set of DUFs assosiated with partical fungi
dufOtherEucaryota = dict()
dufAllSpecies = dict()
ascomycotaDUFs = dict()
basidiomycotaDUFs = dict()
###      declaring empty sets   #####################################################################

dufs = set() #set containing all pfam ids that are DUF
duf_fungi = set()#set containing duf ids thallFungi_setallFungi_setallFungi_setat present in fungithis set is equal dufs set as we are testing dufs in FG fungi
dufFG = set()#set containing duf ids that present in FG (this should be all DUF ids)

duf_planta = set() #set containing duf ids that present also in planta
duf_archaea = set() #set containing duf ids that present also in archaea	
duf_other_eukaryota = set()#set containing duf ids that present also in other eukaryota
duf_animal = set()#set containing duf ids that present also  in animals
duf_bacteria = set()#set containing duf ids that present also in bacteria
duf_unclass = set()#set containing duf ids that present also in unclassified

all_fungi = set()#set containing pfam ids that present in fungi
FGpfam = set()#set containing pfam ids that present in FG - this should be all pfam ids
all_planta = set() #set containing pfam ids that present in planta
all_archaea = set() #set containing pfam ids that present in archaea	
all_other_eukaryota = set()#set containing pfam ids that present inother eukaryota
all_animal = set()#set containing pfam ids that present in animals
all_bacteria = set()#set containing pfam ids that present in bacteria
all_unclass = set()#set containing pfam ids that present in unclassified
allFungi_set = set()
ascomycotaFungi = set()
basidiomycotaFungi = set()
otherFungi = set()
fungiWithDuf_set = set()

pfamInFungi = dict()
allfungiDict = dict()
allAscomycota = dict()
allAscomycotaFungi = set()
allBasidiomycota = dict()
allBasidiomycotaFungi = set()
allOtherFungi = set()

###########################################################################
##########              MAIN PROGRAM             ##########################
###########################################################################
creatingPfamDufsMap_file(infile_pfam27, outfile1, pfam27_to_Duf)	


for line in fh1.readlines():
	line =  line.split('\t')
	pfamAcc = line[0]
	superkingdom = line[2]
	kingdom = line[3]
	phylum = line[4]
	species = line[9].strip()
	
####   Checking all PFAM IDs that are DUFs for taxonomy      ######
	if pfamAcc in pfam27_to_Duf:
		pfamAcc = pfam27_to_Duf[pfamAcc]
		#dufFungiDic.setdefault(pfamAcc,set())#the dictonary where keys are DUFids and value is a set of fungi species containing those DUF ids		
		if species != '':
			dufAllSpecies.setdefault(pfamAcc, set())#creating dict with key DUF id and value is a set of species having this DUF
			if species not in dufAllSpecies[pfamAcc]:
				dufAllSpecies[pfamAcc].add(species) 
				dufOtherEucaryota.setdefault(pfamAcc,set())
		#if pfamAcc == 'DUF3759':
			#fhOut_4.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(pfamAcc, line[1], superkingdom, kingdom, line[4], line[5], line[6], line[7], line[8], species))
		#if pfamAcc == 'DUF2841':
			#fhOut_4.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(pfamAcc, line[1], superkingdom, kingdom, line[4], line[5], line[6], line[7], line[8], species))	
		#if pfamAcc == 'DUF4449':
			#fhOut_4.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(pfamAcc, line[1], superkingdom, kingdom, line[4], line[5], line[6], line[7], line[8], species))
		
		#if pfamAcc not in dufs:
			#dufs.add(pfamAcc)
		if line [2].startswith("Archaea"):
			if pfamAcc not in duf_archaea:
				duf_archaea.add(pfamAcc)
		if line[2].startswith("Eukaryota"):
			if line [3].startswith("Fungi"):
				if species != '':
					fungiWithDuf.setdefault(species, set())
					dufFungiDic.setdefault(pfamAcc,set())#the dictonary where keys are DUFids and value is a set of fungi species containing those DUF ids		
					if species not in dufFungiDic[pfamAcc]:
						dufFungiDic[pfamAcc].add(species)#keys are DUF and value are set of Fungal species
						if line[4].startswith('Ascomycota'):
							ascomycotaDUFs.setdefault(pfamAcc, set())
							if species not in ascomycotaDUFs[pfamAcc]:
								ascomycotaDUFs[pfamAcc].add(species)	
							if species not in ascomycotaFungi:
								ascomycotaFungi.add(species)
						elif line[4].startswith('Basidiomycota'):
							basidiomycotaDUFs.setdefault(pfamAcc, set())
							if species not in basidiomycotaDUFs[pfamAcc]:
								basidiomycotaDUFs[pfamAcc].add(species)
							if species not in basidiomycotaFungi:
								basidiomycotaFungi.add(species)
						else:
							otherFungi.add(species)
					if pfamAcc not in fungiWithDuf[species]:
						fungiWithDuf[species].add(pfamAcc)
					if species not in fungiWithDuf_set:
						fungiWithDuf_set.add(species)
				
					if pfamAcc not in duf_fungi:
						duf_fungi.add(pfamAcc)
			if line [3].startswith("Viridiplantae"):
				if pfamAcc not in duf_planta:
					duf_planta.add(pfamAcc)
				#if pfamAcc == 'DUF3759':
					fhOut_3.write("%s\t%s\n" %(pfamAcc, species))
				#if pfamAcc == 'DUF2841':
					fhOut_3.write("%s\t%s\n" %(pfamAcc, species))	
				#if pfamAcc == 'DUF4449':
					fhOut_3.write("%s\t%s\n" %(pfamAcc, species))	
					
			if line[3].startswith("Metazoa"):
				if pfamAcc not in duf_animal:
					duf_animal.add(pfamAcc)
			if line[3] == '':
				if species not in dufOtherEucaryota[pfamAcc]:
					dufOtherEucaryota[pfamAcc].add(species)
				if pfamAcc not in duf_other_eukaryota:
					duf_other_eukaryota.add(pfamAcc)	
		if line[2].startswith("Bacteria"):
			if pfamAcc not in duf_bacteria:
				duf_bacteria.add(pfamAcc)
		
		if line[2] == '':
			if pfamAcc not in duf_unclass:
				duf_unclass.add(pfamAcc)
				fhOut2.write("%s\t%s\t%s\n" %(line[0], line[1], pfamAcc))
##########    FOR ALL DOMAINS INCLUDING DUFs    #######################################

	if line [2].startswith("Archaea"):
			if pfamAcc not in all_archaea:
				all_archaea.add(pfamAcc)
	if line[2].startswith("Eukaryota"):
		if line [3].startswith("Fungi"):
			if species != '':
				pfamInFungi.setdefault(pfamAcc,set())
				allfungiDict.setdefault(species, set())
				
				if species not in pfamInFungi[pfamAcc]:
					pfamInFungi[pfamAcc].add(species)
					if line[4].startswith('Ascomycota'):
						allAscomycota.setdefault(pfamAcc, set())
						if species not in allAscomycota[pfamAcc]:
							allAscomycota[pfamAcc].add(species)
						if species not in allAscomycotaFungi:	
							allAscomycotaFungi.add(species)
					elif line[4].startswith('Basidiomycota'):
							allBasidiomycota.setdefault(pfamAcc, set())
							if species not in allBasidiomycota[pfamAcc]:
								allBasidiomycota[pfamAcc].add(species)
							if species not in allBasidiomycotaFungi:
								allBasidiomycotaFungi.add(species)	
					else:
							allOtherFungi.add(species)			
										
						
				
					
				if pfamAcc not in allfungiDict[species]:
					allfungiDict[species].add(pfamAcc)
				
				if species not in allFungi_set:
					allFungi_set.add(species)
				
				if pfamAcc not in all_fungi:
					all_fungi.add(pfamAcc)
				
				
		if line [3].startswith("Viridiplantae"):
			if pfamAcc not in all_planta:
				all_planta.add(pfamAcc)
		if line[3].startswith("Metazoa"):
			if pfamAcc not in all_animal:
				all_animal.add(pfamAcc)
		if line[3] == '':
			if pfamAcc not in all_other_eukaryota:
				all_other_eukaryota.add(pfamAcc)	
	if line[2].startswith("Bacteria"):
		if pfamAcc not in all_bacteria:
			all_bacteria.add(pfamAcc)
		
	if line[2] == '':
		if pfamAcc not in all_unclass:
			all_unclass.add(pfamAcc)
			fhOut3.write("%s\t%s\t%s\n" %(line[0], line[1], pfamAcc))
fhOut2.close()
fhOut3.close()
fhOut_3.close()
fhOut_4.close()


###    Calculation of the number of unique fungi species having DUFs that are present in FG (including FG)
duf_diffSpecies = set()
unique_dufOtherEucaryota = set()
for k, v in dufFungiDic.items():
	for n in v:
		if n not in duf_diffSpecies:
			duf_diffSpecies.add(n)
		
for k, v in dufOtherEucaryota.items():
	for n in v:
		if n not in unique_dufOtherEucaryota:
			unique_dufOtherEucaryota.add(n)			
######################################################################################################
#Calculating the number of unique pfam Ids in the taxonomy file results
#######################################################################################################
numberOfPfams = set()
fh = open (infile1, "r")
for line in fh.readlines():
	line = line.split('\t')
	pfamAcc = line[0]
	if pfamAcc in pfam27_to_Duf:
		pfamAcc = pfam27_to_Duf[pfamAcc]
	if pfamAcc not in numberOfPfams:
		numberOfPfams.add(pfamAcc)
		
print "Number of unique pfam IDs in FG classified into Taxid:", len(numberOfPfams)
print "###############################################################################"			
print "number of FG DUFs domain with taxaId:", len(dufs)
print "number of fungi species having DUF ids as FG including FG", len(dufFungiDic)
print "number of DUF Ids present in Fungi", len(duf_fungi)
print "number of DUF ids present also in planta", len(duf_planta)
print "number of DUF ids present also in animala", len(duf_animal)
print "number of DUF Ids present also in Archaea", len(duf_archaea)
print "number of DUF Ids present also in other Eukaryota", len(duf_other_eukaryota)
print "number of DUF IDs present also in Bacteria", len(duf_bacteria)
print "number of DUF IDS present also in unclassified", len(duf_unclass)
duf_fungiOnly = duf_fungi.difference(duf_planta, duf_animal, duf_bacteria, duf_archaea, duf_other_eukaryota)
print "number of DUFs only present in Fungi", len(duf_fungiOnly)
duf_fungiInEukaryota = duf_fungi.difference(duf_bacteria, duf_archaea)
print "fungi with DUFs but without dufs present in bacteria and archaea:", len(duf_fungiInEukaryota)
duf_plantaInEukaryota = duf_planta.difference(duf_bacteria, duf_archaea)
print "plant with DUFs but without dufs present in bacteria and archaea:", len(duf_plantaInEukaryota)
duf_animalsInEukaryota = duf_animal.difference(duf_bacteria, duf_archaea)
print "animal with DUFs but without dufs present in bacteria and archaea:", len(duf_animalsInEukaryota)
duf_other_eukaryotaInEukaryota = duf_other_eukaryota.difference(duf_bacteria, duf_archaea)
print "other eukaryota with DUFs but without dufs present in bacteria and archaea:", len(duf_other_eukaryotaInEukaryota)
print "Number of DUFs only present in Fungi:", duf_fungiOnly
print"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
print "number of different fungi species having PFAM ids as FG including dufs", len(allFungiDic)
print "number of Pfam Ids present in Fungi", len(all_fungi)
print "number of Pfam ids present in planta", len(all_planta)
print "number of Pfam ids present in animala", len(all_animal)
print "number of Pfam Ids present in Archaea", len(all_archaea)
print "number of Pfam Ids present in other Eukaryota", len(all_other_eukaryota)
print "number of Pfam IDs present in Bacteria", len(all_bacteria)
print "number of Pfam IDS present in unclassified", len(all_unclass)
print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
print "##################################################################################"
print "#############        DUF ONLY INTERSECTION OF 2 SETS      ########################"
print "##################################################################################"
print "No FG DUFs domains the same in Fungi and Planta (intersection):", len(duf_fungi.intersection(duf_planta))
print "No FG DUFs domains the same in Fungi and Animals (intersection):", len(duf_fungi.intersection(duf_animal))
print "No FG DUFs domains the same in Fungi and Bacteria (intersection):", len(duf_fungi.intersection(duf_bacteria))
print "No FG DUFs domains the same in Fungi and Archaea:", len(duf_fungi.intersection(duf_archaea))
print "No FG DUFs domains the same in Fungi and undefined (intersection):", len(duf_fungi.intersection(duf_unclass))
print "No FG DUFs domains the same in Fungi and other Eukaryota (intersection):", len(duf_fungi.intersection(duf_other_eukaryota))
print "###########################################################################################"
print "No FG DUFs domains the same in Animal and Planta(intersection):", len(duf_animal.intersection(duf_planta))
print "No FG DUFs domains the same in Animal and Bacteria:", len(duf_animal.intersection(duf_bacteria))
print "No FG DUFs domains the same in Animal and Archaea:", len(duf_animal.intersection(duf_archaea))
print "No FG DUFs domains the same in Animal and undefined (intersection):", len(duf_animal.intersection(duf_unclass))
print "No FG DUFs domains the same in Animal and other Eukaryota (intersection):", len(duf_animal.intersection(duf_other_eukaryota))
print "############################################################################################"
print "No FG DUFs domains the same in Planta and bacteria (intersection):", len(duf_planta.intersection(duf_bacteria))
print "No FG DUFs domains the same in Planta and Archaea:", len(duf_planta.intersection(duf_archaea))
print "No FG DUFs domains the same in Planta and undefined(intersection):", len(duf_planta.intersection(duf_unclass))
print "No FG DUFs domains the same in Planta and other Eukaryota(intersection):", len(duf_planta.intersection(duf_other_eukaryota))
print "############################################################################"
print "No FG DUFs domains the same in bacteria and other Eukaryota (intersection):", len(duf_bacteria.intersection(duf_other_eukaryota))
print "No FG DUFs domains the same in bacteria and undefined (intersection):", len(duf_bacteria.intersection(duf_unclass))
print "No FG DUFs domains the same in bacteria and Archaea (intersection):", len(duf_bacteria.intersection(duf_archaea))
print "#############################################################################################"
print "No FG DUFs domains the same in Archaea and undefined (intersection):", len(duf_archaea.intersection(duf_unclass))
print "No FG DUFs domains the same in Archaea and other Eukaryota(intersection):", len(duf_archaea.intersection(duf_other_eukaryota))
print "#########################################################################"
print "No FG DUFs domains the same in other Eukaryota and undefined(intersection):", len(duf_other_eukaryota.intersection(duf_unclass))
print "##############################################################################################"
print "difference all pfam ids - fungi pfam Ids:", len(numberOfPfams.difference(all_fungi))
print numberOfPfams.difference(all_fungi)
#######################################################################################
############      CHECKING INTERSECTION OF 3 SETS FOR DUFs      #######################
#######################################################################################
duf_FunPlnAnm = set.intersection(duf_fungi, duf_planta, duf_animal) #set of the same DUF ids for Fungi, Planta and Animal
duf_FunPlnOtherEuk = set.intersection(duf_fungi, duf_planta, duf_other_eukaryota)
duf_FunAnmOtherEuk = set.intersection(duf_fungi, duf_animal, duf_other_eukaryota)
duf_PlnAnmOtherEuk = set.intersection(duf_planta, duf_animal, duf_other_eukaryota)
duf_FunPlnBact = set.intersection(duf_fungi, duf_planta, duf_bacteria)#set of the same DUF ids for Fungi, Planta and Bacteria
duf_FunAnmBact = set.intersection(duf_fungi, duf_animal, duf_bacteria)#set of the same DUF ids for Fungi, Animal and Bacteria
duf_PlnAnmBact = set.intersection(duf_planta, duf_animal, duf_bacteria)#set of the same DUF ids for Planta and Animal and Bacteria
print "Number of FG DUF ids common for Fungi, Planta and Metazoa:", len(duf_FunPlnAnm)
print "Number of FG DUF ids common for Fungi, Planta and Other Eukaryota:", len(duf_FunPlnOtherEuk)
print "Number of FG DUF ids common for Fungi, Metazoa and Other Eukaryota:",len(duf_FunAnmOtherEuk)
print "Number of FG DUF ids common for Plant, Metazoa and Other Eukaryota:", len(duf_PlnAnmOtherEuk)
print "Number of FG DUF ids common for Fungi, Plant and Bacteria:", len(duf_FunPlnBact)
print "Number of FG DUF ids common for Fungi, Animal and Bacteria:", len(duf_FunAnmBact)
print "Number of FG DUF ids common for Plant, Animal and Bacteria:", len(duf_PlnAnmBact)
#print "Number of FG DUF ids common for Fungi, Planta, Metazoa and Other Eukaryota:", len(duf_FunPlnAnmOtherEuk)
###################################################################################################################################
##    CHECKING INTERSECTION OF 4 & 5 SETS FOR DUFs#########################################
############################################################################
duf_FunPlnAnmOtherEuk = set.intersection(duf_planta, duf_animal, duf_fungi, duf_other_eukaryota)
duf_FunPlnAnmBact = set.intersection(duf_fungi, duf_planta, duf_animal, duf_bacteria, duf_unclass)
duf_FunPlnBactOtherEuk = set.intersection(duf_fungi, duf_planta, duf_bacteria, duf_other_eukaryota)
duf_AnmPlnBactOtherEuk = set. intersection(duf_animal, duf_planta, duf_bacteria, duf_other_eukaryota)
duf_FunAnmPlnBactOtherEuk = set.intersection(duf_fungi, duf_animal, duf_planta, duf_bacteria, duf_other_eukaryota)
print "Number of FG DUF ids common for Fungi, Planta, Metazoa and Other Eukaryota:", len(duf_FunPlnAnmOtherEuk)
print "Number of FG DUF ids common for Fungi, Planta, Metazoa and Bacteria:", len(duf_FunPlnAnmBact)
print "Number of FG DUF ids common for Fungi, Planta, Bacteria and other Eukaryota:", len(duf_FunPlnBactOtherEuk)
print "Number of FG DUF ids common for Animal, Planta, Bacteria and other Eukaryota:", len(duf_AnmPlnBactOtherEuk)
print "Number of FG DUF ids common for Fungi, Animal, Planta, Bacteria and other Eukaryota:", len(duf_FunAnmPlnBactOtherEuk)
#########################################################################
#      CHECKING INTERSECTION OF 3 SETS FOR ALL PFAM IDs IN FG  #######################
########################################################################
print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
print "##################################################################################"
print "#############        ALL PFAM INTERSECTION OF 2 SETS      ########################"
print "##################################################################################"
print "No FG PFAM domains the same in Fungi and Planta (intersection):", len(all_fungi.intersection(all_planta))
print "No FG PFAM domains the same in Fungi and Animals (intersection):", len(all_fungi.intersection(all_animal))
print "No FG PFAM domains the same in Fungi and Bacteria (intersection):", len(all_fungi.intersection(all_bacteria))
print "No FG PFAM domains the same in Fungi and Archaea:", len(all_fungi.intersection(all_archaea))
print "No FG PFAM domains the same in Fungi and undefined (intersection):", len(all_fungi.intersection(all_unclass))
print "No FG PFAM domains the same in Fungi and other Eukaryota (intersection):", len(all_fungi.intersection(all_other_eukaryota))
print "###########################################################################################"
print "No FG PFAM domains the same in Animal and Planta(intersection):", len(all_animal.intersection(all_planta))
print "No FG PFAM domains the same in Animal and Bacteria:", len(all_animal.intersection(all_bacteria))
print "No FG PFAM domains the same in Animal and Archaea:", len(all_animal.intersection(all_archaea))
print "No FG PFAM domains the same in Animal and undefined (intersection):", len(all_animal.intersection(all_unclass))
print "No FG PFAM domains the same in Animal and other Eukaryota (intersection):", len(all_animal.intersection(all_other_eukaryota))
print "############################################################################################"
print "No FG PFAM domains the same in Planta and bacteria (intersection):", len(all_planta.intersection(all_bacteria))
print "No FG PFAM domains the same in Planta and Archaea:", len(all_planta.intersection(all_archaea))
print "No FG PFAM domains the same in Planta and undefined(intersection):", len(all_planta.intersection(all_unclass))
print "No FG PFAM domains the same in Planta and other Eukaryota(intersection):", len(all_planta.intersection(all_other_eukaryota))
print "############################################################################"
print "No FG PFAM domains the same in bacteria and other Eukaryota (intersection):", len(all_bacteria.intersection(all_other_eukaryota))
print "No FG PFAM domains the same in bacteria and undefined (intersection):", len(all_bacteria.intersection(all_unclass))
print "No FG PFAM domains the same in bacteria and Archaea (intersection):", len(all_bacteria.intersection(all_archaea))
print "#############################################################################################"
print "No FG PFAM domains the same in Archaea and undefined (intersection):", len(all_archaea.intersection(all_unclass))
print "No FG PFAM domains the same in Archaea and other Eukaryota(intersection):", len(all_archaea.intersection(all_other_eukaryota))
print "#########################################################################"
print "No FG PFAM domains the same in other Eukaryota and undefined(intersection):", len(all_other_eukaryota.intersection(all_unclass))
print "##############################################################################################"
print "difference all pfam ids - fungi pfam Ids:", len(numberOfPfams.difference(all_fungi))
print numberOfPfams.difference(all_fungi)
#######################################################################################
############      CHECKING INTERSECTION OF 3 SETS FOR DUFs      #######################
#######################################################################################
FunPlnAnm = set.intersection(all_fungi, all_planta, all_animal) #set of the same DUF ids for Fungi, Planta and Animal
FunPlnOtherEuk = set.intersection(all_fungi, all_planta, all_other_eukaryota)
FunAnmOtherEuk = set.intersection(all_fungi, all_animal, all_other_eukaryota)
PlnAnmOtherEuk = set.intersection(all_planta, all_animal, all_other_eukaryota)
FunPlnBact = set.intersection(all_fungi, all_planta, all_bacteria)#set of the same DUF ids for Fungi, Planta and Bacteria
FunAnmBact = set.intersection(all_fungi, all_animal, all_bacteria)#set of the same DUF ids for Fungi, Animal and Bacteria
PlnAnmBact = set.intersection(all_planta, all_animal, all_bacteria)#set of the same DUF ids for Planta and Animal and Bacteria
FunOtherEukBact = set.intersection(all_fungi, all_other_eukaryota, all_bacteria)
OtherEukPlnBact = set.intersection(all_other_eukaryota, all_planta, all_bacteria)
OtherEukAnmBact = set.intersection(all_other_eukaryota, all_animal, all_bacteria)
print "Number of FG PFAM ids common for Fungi, Planta and Metazoa:", len(FunPlnAnm)
print "Number of FG PFAM ids common for Fungi, Planta and Other Eukaryota:", len(FunPlnOtherEuk)
print "Number of FG PFAM ids common for Fungi, Metazoa and Other Eukaryota:",len(FunAnmOtherEuk)
print "Number of FG PFAM ids common for Plant, Metazoa and Other Eukaryota:", len(PlnAnmOtherEuk)
print "Number of FG PFAM ids common for Fungi, Plant and Bacteria:", len(FunPlnBact)
print "Number of FG PFAM ids common for Fungi, Animal and Bacteria:", len(FunAnmBact)
print "Number of FG PFAM ids common for Plant, Animal and Bacteria:", len(PlnAnmBact)
print "Number of FG PFAM ids common for Fungi, Other Eukaryota and Bacteria:", len(FunOtherEukBact)
print "Number of FG PFAM ids common for Other Eukaryota, Planta and Bacteria:", len(OtherEukPlnBact)
print "Number of FG PFAM ids common for Other Eukaryota, Animal and Bacteria:", len(OtherEukAnmBact)
#print "Number of FG PFAM ids common for Fungi, Planta, Metazoa and Other Eukaryota:", len(FunPlnAnmOtherEuk)
###################################################################################################################################
##    CHECKING INTERSECTION OF 4 & 5 SETS FOR DUFs#########################################
############################################################################
FunPlnAnmOtherEuk = set.intersection(all_planta, all_animal, all_fungi, all_other_eukaryota)
FunPlnAnmBact = set.intersection(all_fungi, all_planta, all_animal, all_bacteria)
FunPlnBactOtherEuk = set.intersection(all_fungi, all_planta, all_bacteria, all_other_eukaryota)
AnmPlnBactOtherEuk = set. intersection(all_animal, all_planta, all_bacteria, all_other_eukaryota)
FunAnmPlnBactOtherEuk = set.intersection(all_fungi, all_animal, all_planta, all_bacteria, all_other_eukaryota)
FunAnmBactOtherEuk = set.intersection(all_fungi, all_animal, all_bacteria, all_other_eukaryota)
print "Number of FG pfam ids common for Fungi, Planta, Metazoa and Other Eukaryota:", len(FunPlnAnmOtherEuk)
print "Number of FG pfam ids common for Fungi, Planta, Metazoa and Bacteria:", len(FunPlnAnmBact)
print "Number of FG pfam ids common for Fungi, Planta, Bacteria and other Eukaryota:", len(FunPlnBactOtherEuk)
print "Number of FG pfam ids common for Fungi, Animal, Bacteria and other Eukaryota:", len(FunAnmBactOtherEuk)
print "Number of FG pfam ids common for Animal, Planta, Bacteria and other Eukaryota:", len(AnmPlnBactOtherEuk)
print "Number of FG pfam ids common for Fungi, Animal, Planta, Bacteria and other Eukaryota:", len(FunAnmPlnBactOtherEuk)
print "#################################################################################################"

dufFungiOnlyDic = dict()
#print 'DUFs in Fungi only with species:'
count = 0
for i in duf_fungiOnly:
	if i in dufFungiDic:
		count = count +1
		#print count, i, dufFungiDic[i]
		dufFungiOnlyDic.setdefault(i, dufFungiDic[i])
		fhOut4.write('%s\t%s\t%s\n' %(count, i, dufFungiOnlyDic[i]))
		
dufFungiOnlySpecies = set()
		
for i in duf_fungiOnly:
	if i in dufFungiDic:
		for e in dufFungiDic[i]:
			if e not in dufFungiOnlySpecies:
				dufFungiOnlySpecies.add(e)
						
	
	
fungiWithDufs = set()
for e in duf_fungiOnly:
	if e in dufFungiDic:
		l = ", ".join(str(i).replace('\r\n', '') for i in dufFungiDic[e])
			#f#or i in range(1, len(fungiDic[e])-1):
			#if value not in fungi:
					#fungi.add()
		newList = l.split(',')
		#print newList
		for element in newList:
			if element not in fungiWithDufs:
				fungiWithDufs.add(element)
			if element == "":
				fungiWithDufs.remove(element)
#print newList
				
##############################################################################################
### CHECKING DUF IDs IF WHAT TYPE FUNGI THEY PRESENT IN MOSTLY:PATHOGENIC OR NOT PATHOGENIC###
##############################################################################################

def checkTheList(infile,out_list, dic):
	duf_freq= dict()
	fh = open (infile, "r")
	index = 1
	numberOfSpecies = set()
	for line in fh.readlines():
		line =  line.split('\t')
		fungi_name = line[0]
		for k in dufFungiOnlyDic:
			if fungi_name in dufFungiOnlyDic[k]:
				out_list.append(k)
				dic.setdefault(k, set())
				if fungi_name not in dic[k]:
					dic[k].add(fungi_name)
			if fungi_name not in numberOfSpecies:
				numberOfSpecies.add(fungi_name)
	print k, numberOfSpecies
					
			
	for i in out_list:
		duf_freq[i] = duf_freq.get(i,0) +1
		
	for k, v in duf_freq.items():
		print k, v
	for k, v in dic.items():
		print k, len(v), v
		
	fh.close()
	return(duf_freq)
	
plant_pathogens_dufs = list()
animal_pathogen_dufs = list()
fungi_pathogen_duf = list()
plant_other_symb_dufs = list()
no_pathogen_dufs = list()

plant_pathogens = dict()
animal_pathogen = dict()
fungi_pathogen = dict()
plant_other_symb = dict()
no_pathogen = dict()	

print "Plant pathogenic_Fungi_with DUFs:"
checkTheList(infile5, plant_pathogens_dufs, plant_pathogens)
#print "Plant pathogenic_Other_Fungi_with DUFs:"	
#checkTheList(infile4, plant_other_symb_dufs, plant_other_symb)
#print "Fungi_Pathogenic_Fungi_with_DUFs:"	
#checkTheList(infile3, fungi_pathogen_duf, fungi_pathogen)
#print "Animal_pathogenic_Fungi_with_DUFs:"	
#checkTheList(infile6, animal_pathogen_dufs, animal_pathogen)
#print "no_pathogenic_Fungi_with_DUFs:"	
#checkTheList(infile2, no_pathogen_dufs, no_pathogen)


