#!/usr/bin/python
import community
from community import *
import networkx as nx
import re
from scipy.stats.stats import pearsonr
import scipy
import matplotlib.pyplot as plt
from colour import Color
print "starting script.."

workDir1 = "/home/ela/Project/pfam(sinceJune2013)"

infile1 = "%s/input/pfamInPhi-base/PHI_FGSG_List" % (workDir1)
infile2 = "%s/output/bigramsFG/bigram_orderNoMatter_12Nov/withDufs/withoutStartEnd/bothBigrams_piSeparated/pfamDomCount_inFG_withUniquePfamNo.txt" % (workDir1)
infile3 = "%s/input/bigramsNetwork/domain_clustering/louvain_clustering/hetero_bigrams_tabDel.txt" % (workDir1)
outfile = "%s/output/fg_bigramsNetwork/bigramOrderNoMatter/heteroBigramsNetwork/increasedVirulence.txt" % (workDir1)
outfile1 ="%s/output/fg_bigramsNetwork/bigramOrderNoMatter/heteroBigramsNetwork/connectedComponentsInfo.txt" % (workDir1)



fh1 = open (infile1, "r")
fh2 = open (infile2, "r")
fh3 = open (infile3, "r")
fhout = open (outfile, "w")
fh_out1 = open (outfile1, "w")

fg_pfam = dict()
domSet = set()	

for line in fh2.xreadlines():
	fields = line.split()
	fgId= fields[0].strip()
	pfamDoms = fields[3:]

	for dom in pfamDoms:
		fg_pfam.setdefault(fgId,set())
		if dom not in fg_pfam[fgId]:
			fg_pfam[fgId].add(dom)

phiFg_InVr = dict()
phiFg_Lethal = dict()
phiFg_LossPath = dict()
phiFg_RedVir = dict()
phiFg_UnaffPath = dict()
phiFg_MixOut = dict()
phiFg_UnaffVir = dict()

pfamInVr = set()
pfamLetal = set()
pfamLossPath = set()
pfamRedVir = set()
pfamUnaffPath = set()
pfamMixOut = set()
pfamUnaffVir = set()

def getPhenoTyDict(a,d,k,v):
	if phenoType.startswith(a):
		d.setdefault(k,v)
     	return len(d)

def pfamWithPhenotype(d,s,t):
	c =0
	for k,v in fg_pfam.items():
		if d.has_key(k):
			c=c+1
			print t, c, k, v
			for e in v:
				if e not in s:
					s.add(e)
			
			#print t, "pfam domains:"
	return ( t+ " pfam domains:", s)


for line in fh1.xreadlines():
	fields = line.split('\t')
	phiBaseId = fields[0].strip()
	fgId_inPhi = fields[1].strip()
	phenoType = fields[3].strip()
	
	lethal = getPhenoTyDict('Lethal', phiFg_Lethal, fgId_inPhi,phiBaseId)
	incVir = getPhenoTyDict('Increased virulence', phiFg_InVr, fgId_inPhi,phiBaseId)
	LossPath = getPhenoTyDict('Loss of pathogenicity', phiFg_LossPath, fgId_inPhi,phiBaseId)
	RedVir = getPhenoTyDict('Reduced virulence', phiFg_RedVir, fgId_inPhi,phiBaseId)
	UnaffPath = getPhenoTyDict('Unaffected pathogenicity', phiFg_UnaffPath, fgId_inPhi,phiBaseId)
	UnaffVir = getPhenoTyDict('Unaffected virulence', phiFg_UnaffVir, fgId_inPhi,phiBaseId)
	MixOut = getPhenoTyDict('mixed outcome', phiFg_MixOut, fgId_inPhi,phiBaseId)

print "Unique No of Lethal phenotypes in FG- source PHI-base:", lethal
print "Unique No of Increase virulence phenotypes in FG- source PHI-base:", incVir
print "Unique No of Increase virulence phenotypes in FG- source PHI-base:", LossPath
print "Unique No of Reduced virulence phenotypes in FG- source PHI-base:", RedVir
print "Unique No of Unaffected pathogenicity phenotypes in FG- source PHI-base:",UnaffPath
print "Unique No of Unaffected virulence phenotypes in FG- source PHI-base:",UnaffVir
print "Unique No of mixed outcome phenotypes in FG- source PHI-base:", MixOut
	

print pfamWithPhenotype(phiFg_InVr,pfamInVr, "Increased virulence")
print pfamWithPhenotype(phiFg_Lethal,pfamLetal, "Lethal")
Lethal= pfamWithPhenotype(phiFg_InVr,pfamInVr, "Lethal")
fhout.write(' '.join(map(str,Lethal)))
print pfamWithPhenotype(phiFg_LossPath,pfamLossPath, "Loss of pathogenicity")
print pfamWithPhenotype(phiFg_RedVir,pfamRedVir, "Reduced virulence")
print pfamWithPhenotype(phiFg_UnaffPath,pfamUnaffPath, "Unaffected pathogenicity")
print pfamWithPhenotype(phiFg_UnaffVir,pfamUnaffVir, "Unaffected virulence")
print pfamWithPhenotype(phiFg_MixOut,pfamMixOut, "mixed outcome")

fhout.close()
G = nx.Graph()
for line in fh3.xreadlines():
	fields = line.split('\t')
	node1 = fields[0].strip()
	node2 = fields[1].strip()
	G.add_edge(node1,node2, weight = float(fields[4].strip()))

fh3.close()
print "nodes:", G.number_of_nodes()
print "edges:", G.number_of_edges()

network_CC = nx.number_connected_components(G)
print "number of connected components:", network_CC

duf_nodes = dict()
#uni_duf_nodes = dict()

i = 1
for k in nx.connected_component_subgraphs(G):
	
	#print i, "CC has", k.number_of_nodes(),"Nodes and ", k.number_of_edges(), "Edges"
	#fh_out1.write("%s\t%s\t%s\t%s" %(i, k.number_of_nodes(), k.number_of_edges(), len(duf_nodes)))
	i = i+1
number = nx.number_connected_components(G)
for n in range(0,(number)):
	H = nx.connected_component_subgraphs(G)[n]
	for k in H.nodes():
		nodesNo = len(H.nodes())
		if k.startswith('DUF'):
			duf_nodes.setdefault(n, list())
			duf_nodes[n].append(k)
			
		
	if duf_nodes.has_key(n):
		#print n, nodesNo,len(duf_nodes[n])
		fh_out1.write("%s\t%s\t%s\t%s\t%s\n" %(n, nodesNo, H.nodes(), len(duf_nodes[n]), duf_nodes[n]))
	else:
		#print n, nodesNo, 0
		fh_out1.write("%s\t%s\t%s\t%s\n" %(n, nodesNo, H.nodes(), 0))	 
fh_out1.close()
#number = nx.number_connected_components(G)
#node_number = 0
#for n in range(0,(number)):
	#H = nx.connected_component_subgraphs(G)[n]
	#size = H.number_of_nodes()
	#i = n+1
	#print "number of nodes in " + str(i) +" component is:", size
	#print "list of nodes:", H.nodes()
	#node_number = node_number + size


whole_part = community.best_partition(G)
print "Modularity of a partition of the whole network:",modularity(whole_part, G)

H = nx.connected_component_subgraphs(G)[0]
part = community.best_partition(H)
print "Partitions: ", float(len(list(part.keys())))
print "Modularity of a partition of the largest CC:",modularity(part, H)

pfamDUF = set()
duf_InVr = set()
duf_Letal = set()
duf_LossPath = set()
duf_RedVir = set()
duf_UnaffPath = set()
duf_MixOut = set()
duf_UnaffVir = set()

for n in H.nodes():
	if n.startswith('DUF'):
		pfamDUF.add(n)

duf_InVr= pfamDUF.intersection(pfamInVr)
duf_Letal = pfamDUF.intersection(pfamLetal)
duf_LossPath = pfamDUF.intersection(pfamLossPath)
duf_RedVir = pfamDUF.intersection(pfamRedVir)
duf_UnaffPath = pfamDUF.intersection(pfamUnaffPath)
duf_MixOut = pfamDUF.intersection(pfamMixOut)
duf_UnaffVir = pfamDUF.intersection(pfamUnaffVir)

for h in H.nodes():
	print h,(part[h])

#drawing graph

size = float(len(set(part.values())))

pos = nx.graphviz_layout(H, prog = "neato")
labels = nx.draw_networkx_labels(H, pos, labels= None, font_size = 8, font_color = 'k', font_family= 'sans-serif', font_weight = 'bold', alpha=5.0)

count = 0
duf_UnaffPath_RedVir = set()
duf_UnaffPath_RedVir = duf_UnaffPath.intersection(duf_RedVir)


for com  in set(part.values()):
	count = count+1
	
	#list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
	list_nodes = list()
	duf_nodes = list()
	colorList = list()
	shapeList = list()
	
	for nodes in part.keys():
		if part[nodes] == 0:
			list_nodes.append(nodes)
		
			
			if nodes in pfamDUF:
				#duf_nodes.append(nodes)
				if nodes in duf_UnaffPath:
					if nodes in duf_UnaffPath_RedVir:
						colorList.append('saddlebrown')
					else:
						colorList.append('c')
				elif nodes in duf_RedVir:
					colorList.append('deeppink')
				elif nodes in duf_Letal:
					colorList.append('chartreuse')
				else:
					colorList.append('red')
			
			
			else:
				colorList.append(str(com/size))	
				
			
	#color = str(count/size)		
	nx.draw_networkx_nodes(H , pos, list_nodes, with_labels = True, node_size = 150, node_color = colorList, alpha = 1.0)
	#nx.draw_networkx_nodes(H , pos, duf_nodes, with_labels = True, node_size = 150, node_color = colorList, node_shape = 's',alpha = 1.0)
nx.draw_networkx_edges(H, pos, alpha=0.5)
plt.show()
print list_nodes		

