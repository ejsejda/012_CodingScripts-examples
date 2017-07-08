from numpy import *
from numpy.lib.io import loadtxt
import networkx as nx

#setting project directory
projectDir1="/home/ela/Project/PHI-base4.1-analysis/output"

#colour and size function
def getColorandSize(nodes, d, s):
	colorList=list()
    	sizeList=list()
     	for n in nodes:
		if 'Plant Pathogen' in d[n] and 'Animal Pathogen' in d[n]:
            		c='y'
		elif 'Plant Pathogen' in d[n] and 'Other Pathogen' in d[n]:
			c = 'r'
		elif 'Animal Pathogen' in d[n] and  'Other Pathogen' in d[n]:
			c = 'b'

        	elif 'Plant Pathogen' in d[n]:
            		c='g'
        
        	elif 'Animal Pathogen' in d[n]:
            		c='c'
        
        	elif 'lethal' in d[n]:
            		c='w'
        
        	elif 'chem target' in d[n]:
            		c='orange'
        
        	elif 'Other Pathogen' in d[n]:
            		c='m'
        
        	colorList.append(c)
             
        
        	if s.has_key(n):
            		sizeList.append(int(s[n]))
        	else:
            		sizeList.append(0)
            
    	return colorList, sizeList
	
        


fileName1="MCL_1.6_Output.txt" #% (projectDir1)
fileName2="phiBaseFile_withHostTaxaGrouping_noPHIid.txt" #% (projectDir1)
fileName3="pathClassOnGeneId.csv" #% (projectDir1)
#Read the cluster information into an array

data1=loadtxt(fileName1, dtype='S')
data2=loadtxt(fileName2, dtype='S', delimiter=';')
#d=dict([(e.split(';')[0], e.split(';')[7]) for e in data2])
d=dict()
for e in data2:
    k='PHI:'+e[0]
    v=e[7]
    d.setdefault(k, set()).add(v)



data3=loadtxt(fileName3, dtype='S', delimiter=';')
s=dict([('PHI:'+e[0], e[1]) for e in data3])

#Create an empty dictionary
cluster=dict()
phi=dict()
G=nx.Graph()

#For each element of the array, that is for each line of the file

for row in data1: 
    #Get the cluster Id
    idCluster=int(row[0])-1
    #Get the information
    info1=row[1].split('|')[0]
   
    if not cluster.has_key(idCluster):
        cluster[idCluster]=list()
       
    cluster[idCluster].append(info1)

for idCluster in arange(len(cluster)):
	for g1 in cluster[idCluster]:
		for g2 in cluster[idCluster]:
			if g1 != g2:
				G.add_edge(g1,g2)
   
              


colorList, sizeList = getColorandSize(G.nodes(), d, s)
    
import matplotlib.pyplot as plt
plt.figure(1,figsize=(12,12))


# layout graphs with positions using graphviz neato
pos=nx.graphviz_layout(G,prog="neato")
C=nx.connected_component_subgraphs(G)
#display labels, it is used for single cluster visualisation
#nx.draw(C[1], with_labels=True, node_size=array(sizeList)*1800, node_color=colorList, alpha=1.0)
     
#display all clusters for network G without labels
nx.draw(G, pos,with_labels=False, node_size=array(sizeList)*50,node_color=colorList, alpha=1.0)

plt.show()

     
