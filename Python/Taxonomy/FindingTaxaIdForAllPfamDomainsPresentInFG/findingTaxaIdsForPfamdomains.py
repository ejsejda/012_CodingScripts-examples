#!/usr/bin/python
###############################################################################################
## The aims of this program is to find taxa id for each pfam domain present in 
## Fusariom graminearum proteome using information from UniProt taxonomy via SPARQL query #####
###############################################################################################
import sys
import traceback
import urllib2
from urllib2 import urlopen as urlopener 
from SPARQLWrapper import SPARQLWrapper, JSON
import time

### In this setting all input files: infile1 and infile2 needs to be in the same directory as the python code, you can adjust it by providing the pathway to those files    ######
infile1 = "overlappingSolvedGoldfile.txt" #raw file of pfam domain in FG proteins in order
infile2 = "allPfamDomains_inFG_VS_pfamDomains_1st_run.txt"
outfile1 = "results_findingTaxaIdForPfamDomains.txt" 
fh1 = open (infile1, "r") 
fh2 = open (outfile1, "w")

######################################################################################
# DECLARATION OF GLOBAL DICTIONARIES AND SETS
######################################################################################
pfam_to_Duf = dict()
duf_count = set()#set of all DUFs in FG proteome and the whole set includes 314 DUFs
pfamId_Desc = dict() #dictionary where keys are PFAM or DUF (where present) domains and values are the description of the domain
fg_prot = dict()#dictionary where keys are FGSG_Id and value is a list of total pfam domain within given FG protein, there are 8478 	FG proteins that have PFAM(including DUF  where exist) domains
fg_uniqPfam = dict()#dictionary where keys are FGSG_Id and value is a list of unique pfam domain within given FG protein
pfamIds = dict()#dictionary where keys are pfam/or DUF ids and the value are total number of FGSG_Ids where the given PFAM/or DUF is 	present, so value can be several the same FGSG_ids
OneDufDomFG= set()
uniqDufs_InOneDufDomFG = set()
pfamNum = dict()

allPfams = set()
testPfams = set()

#######################################################################################
q1="""
		PREFIX up:<http://purl.uniprot.org/core/> 
		PREFIX taxon:<http://purl.uniprot.org/taxonomy/> 
		PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
		PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
		SELECT DISTINCT ?taxon ?name 
		WHERE
		{
		?protein a up:Protein . 
		?protein up:organism ?taxon .
		?protein rdfs:seeAlso <http://purl.uniprot.org/pfam/&&&> .
		?taxon a up:Taxon .
		?taxon up:scientificName ?name .
	

		} """

q2 = """
		PREFIX up:<http://purl.uniprot.org/core/>
		PREFIX taxon:<http://purl.uniprot.org/taxonomy/>
		SELECT ?name ?sub ?names ?rank 
 		WHERE 
		{ 
 		<http://purl.uniprot.org/taxonomy/&&&> up:scientificName ?name .
		?sub up:rank ?rank .
		<http://purl.uniprot.org/taxonomy/&&&> rdfs:subClassOf+ ?sub .
		?sub up:scientificName ?names .
		}"""
###############################################################################################
class SPARQLWrapper1(SPARQLWrapper):
	def _query(self):
		if self.timeout:
			socket.setdefaulttimeout(self.timeout)
		request = self._createRequest()
		tmp = request.get_full_url().replace("format=json", "format=srj")
		request = urllib2.Request(tmp)
		try:
			response = urlopener(request)
			return response, self.returnFormat
		except urllib2.HTTPError, e:
			if e.code == 400:
				raise QueryBadFormed(e.read())
			elif e.code == 404:
				raise EndPointNotFound(e.read())
			elif e.code == 500:
				raise EndPointInternalError(e.read())
			else:
				print traceback.format_exc()
				raise e
###################################################################################################
	
def runSparql(query, variable):
	sparql = SPARQLWrapper1("http://beta.sparql.uniprot.org")
	s = query.replace('&&&', variable)
	#print s
	#sys.exit()
	a = True
	repeat_count = 0
	while a:
		try:
			sparql.setQuery(s)
			sparql.setReturnFormat('json')
			results = sparql.query().convert()
		except KeyboardInterrupt:
			raise	
		except:
			time.sleep(3)
			repeat_count+=1
			if  repeat_count > 20:
				a = False
				print "Failed_completely"
			continue
		a = False	
	return results

#####################################################################################################
def readingGoldFile():

	for line in fh1.readlines():
		line = line.split('\t')
		fgId = line[0].strip()#first column in the file1 is FGSG_id
		pfamAcc = line[1][:7].strip()#second column in the file1 is pfam accesion number and here we remove everything after 7 character of pfam accesion no
		pfamDesc = line[2].strip()
		pfamStart = line[3][:-2].strip()
		pfamEnd = line[4][:-3].strip()
		pfamInfo = str(pfamAcc) +"_" + str(pfamStart) + "-" + str(pfamEnd)#concatenating start and end domain into pfamAccesion
		pfamId_Desc.setdefault(pfamAcc, pfamDesc)
		fg_prot.setdefault(fgId,list())
		fg_uniqPfam.setdefault(fgId,set())
		fg_prot[fgId].append(pfamAcc)
		fg_uniqPfam[fgId].add(pfamAcc)
		pfamIds.setdefault(pfamAcc,list())
		pfamIds[pfamAcc].append(fgId)
	return(pfamId_Desc, fg_prot, fg_uniqPfam, pfamIds, duf_count)
#######################################################################################################

def compare_Two_colums(file):
	fh = open (file, "r")
	
	for line in fh.xreadlines():
		#if line.isspace():
			#pass
		#else:
		line = line.split('\t')
		allPfamID = line[0].strip()
		testPfamID = line[1].strip()
		#print allPfamID, testPfamID
		if allPfamID not in allPfams:
			allPfams.add(allPfamID)
		if line[1].strip():#needs to be as two columns are different size, and column 2 is shorter
			if testPfamID not in testPfams:
				testPfams.add(testPfamID)
	return(allPfams, testPfams)
##########################################################################################################
#####              MAIN PROGRAM      #####################################################################
##########################################################################################################
readingGoldFile()
domains = list()
compare_Two_colums(infile2)
needsTesting = set()
print "Number of PFAM in FG", len(allPfams)
print "Number of PFAM tested for taxid so far", len(testPfams)
needsTesting = allPfams.difference(testPfams)
print 'not assigned taxid:', len(needsTesting)

for k in pfamIds.items():
	if k not in domains:
		domains.append(k)
print len(domains)
d = 0

for j in range(0, len(domains)-1):
	pf = domains[j]
	results = runSparql(q1, str(pf[0]))
		
	for i in range(0,len(pf)-1):
		d = i
		web1 = "http://purl.uniprot.org/uniprot/"
		web2 = "http://purl.uniprot.org/taxonomy/"
		
		for result in results["results"]["bindings"]:
			taxid = result["taxon"]['value'].replace(web2, "")
			results1 = runSparql(q2, taxid)
			kingdom = ""
			superkingdom = ""
			phylum = ""
			subclass = ""
			order = ""
			family =""
			genus = ""
			species = ""
			for result1 in results1["results"]["bindings"]:
				rank = result1["rank"]['value']
				if "Kingdom" in rank:
					kingdom = result1["names"]['value']
				if "Superkingdom" in rank:
					superkingdom = result1["names"]['value']
				
				if "Phylum" in rank:
					phylum = result1["names"]['value']
				
				if "Subclass" in rank:
					subclass = result1["names"]['value']
				if "Order" in rank:
					order = result1["names"]['value']
				
				if "Family" in rank:
					family = result1["names"]['value']
				if "Genus" in rank:
					genus = result1["names"]['value']
				if "Species" in rank:
					species = result1["names"]['value']
						
			print pf[i], result["taxon"]['value'].replace(web2, ""),superkingdom, kingdom, phylum, species
			fh2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(pf[i],result["taxon"]['value'].replace(web2, ""), superkingdom, kingdom, phylum, subclass, order, family, genus, species))

