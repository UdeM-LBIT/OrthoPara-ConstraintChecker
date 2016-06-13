 #!/usr/bin/python
 
import httplib2, sys
import os
import random
import json
import time
 
http = httplib2.Http(".cache")
 
genetreeids = set()
maxNbTrees = 5000
probTree = 0.2

idsfile = open('/u/lafonman/data/ensembl/gene_tree_ids.txt')
for line in idsfile:
  
  if random.random() <= probTree and len(genetreeids) < maxNbTrees:
  
    pz = line.split('\t')
    
    if len(pz) > 0:
      genetreeids.add(pz[0])
idsfile.close()


seqsByTaxon = {}

def visit_json_elem(jsonElem):

	if 'children' in jsonElem:
		for c in jsonElem['children']:
			visit_json_elem(c)
	else:
		tax = jsonElem['taxonomy']['scientific_name']
		if tax not in seqsByTaxon:
			seqsByTaxon[tax] = []
		obj = {}
		obj['id'] = jsonElem['id']['accession']
		obj['seq'] = jsonElem['sequence']['mol_seq']['seq']
		seqsByTaxon[tax].append( obj )


for treeid in genetreeids:
  
  if os.path.isfile(treeid + "-strict.blast-graph"):
    print treeid + "-strict.blast-graph EXISTS"
    continue
  
  print "TREEID=" + treeid
  #treeid = "ENSGT00390000013823"
  server = "http://beta.rest.ensembl.org"
  ext = "/genetree/id/" + treeid + "?"
  resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
  
  if not resp.status == 200:
	  print "Invalid response: ", resp.status
	  continue	#evil continue


  seqsByTaxon.clear()
  seqsByTaxon = {}

  decoded = json.loads(content)


  tree = decoded['tree']

  visit_json_elem(tree)
  
  if len(seqsByTaxon) < 10:
    print treeid + " has only " + str(len(seqsByTaxon)) + " taxa, skipping"
    continue

  command = "proteinortho5.pl "

  for tax in seqsByTaxon:
	  fname = '/u/lafonman/data/treeseqs/' + tax.replace(' ', '_') + '.fasta'

	  command += '"' + fname + '" '
	  f = open(fname, 'w')

	  for obj in seqsByTaxon[tax]:
		  f.write('>' + obj['id'] + '\n' + obj['seq'] + '\n')

	  f.close()


  command += " -graph -clean "

  #do loose to strict extraction with proteinortho
  commands = {}
  commands['ultraloose'] = command + " -identity=10 -cov=10 -conn=0.001 -sim=0.1"
  commands['loose'] = command + " -identity=10 -cov=25 -conn=0.001 -sim=0.25"
  commands['default'] = command
  commands['strict'] = command + " -identity=50 -cov=75 -conn=0.2 -sim=0.95"
  commands['ultrastrict'] = command + " -identity=50 -cov=75 -conn=0.3 -sim=0.99"


  for c in commands:
    curcommand = commands[c] + " -project=" + treeid + "-" + c
    print "EXECUTING " + curcommand
    os.system(curcommand)
    
  time.sleep(1)


