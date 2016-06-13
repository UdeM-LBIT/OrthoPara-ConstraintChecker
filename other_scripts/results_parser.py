 #!/usr/bin/python
 
import httplib2, sys
import os
import subprocess
import random
import json
import time
 
 


indir = "/u/lafonman/data/constraint_checker_results/"
indirEnsemblRels = indir + "ensembl_rels/"


suffixes = ('ultraloose', 'loose', 'default', 'strict', 'ultrastrict')

nbFiles = 0

statsBySingleSuffix = {}

for suffix in suffixes:
  statsBySingleSuffix[suffix] = {"NBISSAT" : 0, "NBISCONS" : 0}
statsBySingleSuffix["NBHASSAT"] = 0
statsBySingleSuffix["NBHASCONS"] = 0

statsByDoubleSuffix = {}
#statsByDoubleSuffix["NBHASSAT"] = 0
#statsByDoubleSuffix["NBHASCONS"] = 0

#-------------------------------------------------------------------------
def getRelsFromString(relstr):
  retrels = set()
  
  relz = relstr.split(',')
      
  for relx in relz:
    if relx != '':	
      rel = relx.split(':')
      #print relx
      #print rel
      retrels.add( (rel[0], rel[1]) )
  
  return retrels
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def getEnsemblRelationships(relfile):
  f = open(relfile)
  
  ret = {}
  
  for line in f:
    line = line.replace('\n', '')
    
    currels = set()
    if '=' in line:
      pz = line.split('=')
      
      ret[pz[0].upper()] = getRelsFromString(pz[1])
    
  f.close()  
  
  return ret
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
def compareWithEnsemblRelationships(myorthologies, myparalogies, ensemblrels):
  
  badorthos = set()
  badparas = set()
  badorthosindubious = set()
  
  ensorth = ensemblrels["ORTHOLOGS"]
  enspara = ensemblrels["PARALOGS"]
  ensdub = ensemblrels["DUBIOUS"]
  
  for orth in myorthologies:
    if not orth in ensorth:
      badorthos.add(orth)
      if orth in ensdub:
	badorthosindubious.add(orth)
      
  for para in myparalogies:
    if not para in enspara and not para in ensdub:
      badparas.add(para)
      
  
  return (badorthos, badparas, badorthosindubious)
#-------------------------------------------------------------------------




for filename in os.listdir(indir):
    if filename.endswith(".satcons"):
	nbFiles+=1
	curDoubleSuffix = ''
        infile = open(indir + filename)
        
        ensemblrels = getEnsemblRelationships(indirEnsemblRels + filename.replace(".satcons", "") + ".ensrels")
        
        hasSingleSAT = False
        hasSingleCONS = False
        curModeNbGenes = 0
        curModeNbOrthologs = 0
        curModeNbParalogs = 0
        curModeOrthoStr = ''
        curModeParaStr = ''
        isCurModeSAT = False
        isCurModeCONS = False
        
        for line in infile:
	  line = line.replace('\n', '')
	  
	  
	  if curDoubleSuffix == '':
	    if line.startswith("<") and not line.startswith("</"):
	      curDoubleSuffix = line.replace("<", "").replace(">", "")
	      
	      if not curDoubleSuffix in statsByDoubleSuffix:
		statsByDoubleSuffix[curDoubleSuffix] = {}
	      
	    elif line.startswith("ISSAT"):
	      pz = line.split("=")
	      issat = pz[1]
	      sz = pz[0].split("-")
	      suffix = sz[1]
	      
	      if issat == "1":
		statsBySingleSuffix[suffix]["NBISSAT"] += 1
		hasSingleSAT = True
	    elif line.startswith("ISCONS"):
	      pz = line.split("=")
	      iscons = pz[1]
	      sz = pz[0].split("-")
	      suffix = sz[1]
	      
	      if iscons == "1":
		statsBySingleSuffix[suffix]["NBISCONS"] += 1
		hasSingleCONS = True
	  #if curDoubleSuffix is set, we're in a mixup of the kind loose-strict
	  else:		
	    if line.startswith("</"):
	      
	      if isCurModeSAT and curModeNbGenes > 0:
		ttlRels = float(curModeNbGenes * (curModeNbGenes - 1)) / 2.0
		pctRels = float(curModeNbOrthologs + curModeNbParalogs)/ttlRels
		
		#if pctRels > 1.0:
		#  print filename
		#  print "satrels=", curModeNbGenes, " choose 2 = ", ttlRels, " nbrels=", curModeNbOrthologs + curModeNbParalogs, " pct=", pctRels
		if not "SUM_PCT_REL_SAT" in statsByDoubleSuffix[curDoubleSuffix]:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_REL_SAT"] = 0.0
		  
		statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_REL_SAT"] += pctRels
		statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_REL_SAT"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_REL_SAT"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISSAT"])
		
		
		#count how many times single has no SAT but combined has SAT
		if not hasSingleSAT:
		  if not "NB_SAT_WHEN_NO_SINGLE_SAT" in statsByDoubleSuffix[curDoubleSuffix]:
		    statsByDoubleSuffix[curDoubleSuffix]["NB_SAT_WHEN_NO_SINGLE_SAT"] = 0
		  statsByDoubleSuffix[curDoubleSuffix]["NB_SAT_WHEN_NO_SINGLE_SAT"] += 1
		
		myorthologies = getRelsFromString(curModeOrthoStr)
		myparalogies = getRelsFromString(curModeParaStr)
		vsensembl = compareWithEnsemblRelationships(myorthologies, myparalogies, ensemblrels)
		
		#print len(vsensembl[0]), "/", len(myorthologies), "  W/O DUBIOUS : ", len(vsensembl[0]) - len(vsensembl[2]), "/", len(myorthologies), " P: ", len(vsensembl[1]), "/", len(myparalogies), " BAD ORTHOS IN DUBIOUS=", len(vsensembl[2])
		
		#compute the average percentage of orthologies/paralogies that don't match those of Ensembl
		if not "SUM_PCT_BADORTH_SAT" in statsByDoubleSuffix[curDoubleSuffix]:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_SAT"] = 0.0
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADPARA_SAT"] = 0.0
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_DUB_SAT"] = 0.0
		
		
		
		
		if len(myorthologies) > 0:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_SAT"] += float(len(vsensembl[0])) / float(len(myorthologies))
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_DUB_SAT"] += float(len(vsensembl[0]) - len(vsensembl[2])) / float(len(myorthologies))
		  statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_BADORTH_SAT"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_SAT"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISSAT"])
		  statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_BADORTH_DUB_SAT"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_DUB_SAT"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISSAT"])
		if len(myparalogies) > 0:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADPARA_SAT"] += float(len(vsensembl[1])) / float(len(myparalogies))
		  statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_BADPARA_SAT"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADPARA_SAT"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISSAT"])
		
		
	      if isCurModeCONS and curModeNbGenes > 0:
		ttlRels = float(curModeNbGenes * (curModeNbGenes - 1)) / 2.0
		pctRels = float(curModeNbOrthologs + curModeNbParalogs)/ttlRels
		
		#print "consrels=", pctRels
		if not "SUM_PCT_REL_CONS" in statsByDoubleSuffix[curDoubleSuffix]:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_REL_CONS"] = 0.0
		statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_REL_CONS"] += pctRels
		statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_REL_CONS"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_REL_CONS"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISCONS"])
	      
		#count how many times single has no CONS but combined has CONS
		if not hasSingleCONS:
		  if not "NB_CONS_WHEN_NO_SINGLE_CONS" in statsByDoubleSuffix[curDoubleSuffix]:
		    statsByDoubleSuffix[curDoubleSuffix]["NB_CONS_WHEN_NO_SINGLE_CONS"] = 0
		  statsByDoubleSuffix[curDoubleSuffix]["NB_CONS_WHEN_NO_SINGLE_CONS"] += 1
	      
		
		myorthologies = getRelsFromString(curModeOrthoStr)
		myparalogies = getRelsFromString(curModeParaStr)
		vsensembl = compareWithEnsemblRelationships(myorthologies, myparalogies, ensemblrels)
		#compute the average percentage of orthologies/paralogies that don't match those of Ensembl
		if not "SUM_PCT_BADORTH_CONS" in statsByDoubleSuffix[curDoubleSuffix]:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_CONS"] = 0.0
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADPARA_CONS"] = 0.0
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_DUB_CONS"] = 0.0
		if len(myorthologies) > 0:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_CONS"] += float(len(vsensembl[0])) / float(len(myorthologies))
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_DUB_CONS"] += float(len(vsensembl[0]) - len(vsensembl[2])) / float(len(myorthologies))
		  statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_BADORTH_CONS"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_CONS"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISCONS"])
		  statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_BADORTH_DUB_CONS"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADORTH_DUB_CONS"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISCONS"])
		if len(myparalogies) > 0:
		  statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADPARA_CONS"] += float(len(vsensembl[1])) / float(len(myparalogies))
		  statsByDoubleSuffix[curDoubleSuffix]["AVG_PCT_BADPARA_CONS"] = statsByDoubleSuffix[curDoubleSuffix]["SUM_PCT_BADPARA_CONS"]/float(statsByDoubleSuffix[curDoubleSuffix]["NBISCONS"])
		
		
		
		
	      
	      curDoubleSuffix = ''
	      curModeNbGenes = 0
	      curModeNbOrthologs = 0
	      curModeNbParalogs = 0
	      isCurModeSAT = False
	      isCurModeCONS = False
	      curModeOrthoStr = ''
	      curModeParaStr = ''
	      
	      
	    elif line.startswith("ISSAT") and line.endswith("1"):
	      
		if not "NBISSAT" in statsByDoubleSuffix[curDoubleSuffix]:
		  statsByDoubleSuffix[curDoubleSuffix]["NBISSAT"] = 0
		statsByDoubleSuffix[curDoubleSuffix]["NBISSAT"] += 1
		isCurModeSAT = True
		
	    elif line.startswith("ISCONS") and line.endswith("1"):
	      
	      if not "NBISCONS" in statsByDoubleSuffix[curDoubleSuffix]:
		statsByDoubleSuffix[curDoubleSuffix]["NBISCONS"] = 0
	      statsByDoubleSuffix[curDoubleSuffix]["NBISCONS"] += 1
	      
	      isCurModeCONS = True
	      
	    elif line.startswith("NBGENES="):
	      curModeNbGenes = int(line.split("=")[1])
		
	    elif line.startswith("NBORTHOLOGS="):
	      curModeNbOrthologs = int(line.split("=")[1])
	      
	    elif line.startswith("ORTHOLOGS="):
	      curModeOrthoStr = line.split("=")[1]
	    
	    elif line.startswith("PARALOGS="):
	      curModeParaStr = line.split("=")[1]
	    
	    elif line.startswith("NBPARALOGS="):
	      curModeNbParalogs = int(line.split("=")[1])
	      
	    
	      
	infile.close()
	if hasSingleSAT:
	  statsBySingleSuffix["NBHASSAT"] += 1
	if hasSingleCONS:
	  statsBySingleSuffix["NBHASCONS"] += 1
	
	
	
	
print "NBTREES=" + str(nbFiles)

for suffix in statsBySingleSuffix:
  print "--------" + suffix + "-----------"
  
  if isinstance(statsBySingleSuffix[suffix], int):
    print suffix, statsBySingleSuffix[suffix]
  else:
    for key in statsBySingleSuffix[suffix]:
      print key + "=" + str(statsBySingleSuffix[suffix][key])

  print "\n\n"

for suffix in sorted(statsByDoubleSuffix.keys()):
  print "--------" + suffix + "-----------"
  
  if isinstance(statsByDoubleSuffix[suffix], int):
    print suffix, statsByDoubleSuffix[suffix]
  else:
    for key in sorted(statsByDoubleSuffix[suffix].keys()):
      print key + "=" + str(statsByDoubleSuffix[suffix][key])

  print "\n\n"
