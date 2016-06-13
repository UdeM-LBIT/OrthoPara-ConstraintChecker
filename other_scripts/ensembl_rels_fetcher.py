 #!/usr/bin/python
 
import httplib2, sys
import os
import random
import json
import time

from ete2 import TreeNode
 
http = httplib2.Http(".cache")
 
genetreeids = set()
maxNbTrees = 5000
probTree = 0.2

indir = '/u/lafonman/data/constraint_checker_results/'
outdir = '/u/lafonman/data/constraint_checker_results/ensembl_rels/'



def visit_json_elem(jsonElem, node, genesNodeMappingToFill):

  if 'events' in jsonElem:
    node.add_feature('event', jsonElem['events']['type'])
    node.name = jsonElem['events']['type']
  if 'children' in jsonElem:
    for c in jsonElem['children']:
      child = node.add_child()
      visit_json_elem(c, child, genesNodeMappingToFill)
  else:
    node.name = jsonElem['id']['accession']
    genesNodeMappingToFill[node.name] = node



for filename in os.listdir(indir):
    if filename.endswith(".satcons"):
      treeid = filename.replace(".satcons", "")
  
      outfilename = outdir + treeid + ".ensrels"
  
      if os.path.isfile(outfilename):
	print outfilename + " exists, skipping"
	continue
      
      print "TREEID=" + treeid
      #treeid = "ENSGT00390000013823"
      server = "http://beta.rest.ensembl.org"
      ext = "/genetree/id/" + treeid + "?"
      resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
      
      if not resp.status == 200:
	print "Invalid response: ", resp.status
	continue	#evil continue


      decoded = json.loads(content)

      tree = decoded['tree']

      geneNodeMap = {}
      root = TreeNode()
      visit_json_elem(tree, root, geneNodeMap)
      
      orthologs = set()
      paralogs = set()
      paralogs_dubious = set()
      
      for g1 in geneNodeMap:
	for g2 in geneNodeMap:
	  if g1 != g2:
	    n1 = geneNodeMap[g1]
	    n2 = geneNodeMap[g2]
	    
	    lca = n1.get_common_ancestor(n2)
	    
	    if lca.name == 'duplication' or lca.name == 'gene_split':
	      paralogs.add( (g1, g2) )
	    elif lca.name == 'dubious':
	      paralogs_dubious.add( (g1, g2) )
	    elif lca.name == 'speciation':
	      orthologs.add( (g1, g2) )
	    else:
	      print "UNDEFINED REL ", g1, g2
	      print "LNAME=", lca.name
	      print lca.features
	      print root
	      sys.exit()
      
      outfile = open(outfilename, 'w')
      
      orthostr = ""
      for ort in orthologs:
	if orthostr != '':
	  orthostr += ","
	orthostr += ort[0] + ":" + ort[1]
	
      parastr = ""
      for par in paralogs:
	if parastr != '':
	  parastr += ","
	parastr += par[0] + ":" + par[1]
	
      dubstr = ""
      for dub in paralogs_dubious:
	if dubstr != '':
	  dubstr += ","
	dubstr += dub[0] + ":" + dub[1]
      
      outfile.write("ORTHOLOGS=" + orthostr + '\n')
      outfile.write("PARALOGS=" + parastr + '\n')
      outfile.write("DUBIOUS=" + dubstr + '\n')
      
      outfile.close()
      
      time.sleep(1)
      


