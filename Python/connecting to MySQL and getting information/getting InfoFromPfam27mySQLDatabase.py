#!/usr/bin/python
# PROGRAM FOR CONNECTION TO MySQL DATABASE pfamdb created on babvs38 WHERE WE STORE TABLES FROM PFAM 27.0. HERE WE ARE INTERESTED IN TWO TABLES: pfamA_ncbi and ncbi_taxonomy AND WE JOIN THE TABLES ON THE SAME ncbi_taxid. THEN RESULTS ARE WRITEN INTO THE FILE -file1

import MySQLdb
import sys

workDir1 = "/home/ela/Project/pfam(sinceJune2013)"
file1 = "%s/output/pfam27_taxonomy/mySql_query_sep2014.txt" % (workDir1)

fh1 = open(file1, "w")

#create the connection to pfamdb on babvs38 and checking the connection
try:
	myPfamDb = MySQLdb.connect(host = "babvs38.rothamsted.ac.uk", db="pfam27_db", user= "ela", passwd="*****" )#password need to be inputed if you would like to re-run code
	print " successful conection"
except StandardError, e:
	print "can't connect to the database", e
	sys.exit()
#create a cursor
cursor = myPfamDb.cursor()
cursor.execute("select pfamA_acc, species, taxonomy from pfamA_ncbi left join ncbi_taxonomy on pfamA_ncbi.ncbi_taxid = ncbi_taxonomy.ncbi_taxid")
rows = cursor.fetchall()
count = 0
for columns in rows:
	l = list(columns)
	count = count +1
	#print l
	fh1.write("%s\t%s\n" %(count,l))
fh1.close()
print count