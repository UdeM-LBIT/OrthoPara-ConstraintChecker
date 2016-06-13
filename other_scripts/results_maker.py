 #!/usr/bin/python
 
import httplib2, sys
import os
import subprocess
import random
import json
import time
 
http = httplib2.Http(".cache")
 
genetreeids = set()

idsfile = open('/u/lafonman/data/ensembl/gene_tree_ids.txt')
for line in idsfile:
  
  pz = line.split('\t')
    
  if len(pz) > 0:
    genetreeids.add(pz[0])
idsfile.close()

#yup, it's all hardcoded
#this is taken from 
# http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-compara/scripts/pipeline/species_tree_blength.nh?root=ensembl&view=co
speciesTreeStr = '((((((((((((((((((((((Pan_troglodytes:0.006667,Homo_sapiens:0.0067):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,Macaca_mulatta:0.037471):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Ictidomys_tridecemlineatus:0.225629):0.01015,Cavia_porcellus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,(((((((Ailuropoda_melanoleuca:0.025614,Mustela_putorius_furo:0.0256):0.0256145,Canis_familiaris:0.051229):0.051229,Felis_catus:0.098612):0.049845,Equus_caballus:0.109397):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508,(((Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,Ovis_aries:0.061796):0.061796):0.025153,Sus_scrofa:0.107275):0.0201675,Vicugna_pacos:0.079):0.0201675):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,((Macropus_eugenii:0.101004,Sarcophilus_harrisii:0.101004):0.021004,Monodelphis_domestica:0.125686):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,(((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,(Ficedula_albicollis:0.085771,Taeniopygia_guttata:0.085771):0.085771):0.199223,Pelodiscus_sinensis:0.489241):0.105143,Anolis_carolinensis:0.4989):0.17):0.149,Xenopus_tropicalis:0.855573):0.155677,Latimeria_chalumnae:0.155677):0.550693,(((((((Xiphophorus_maculatus:0.1204925,Oryzias_latipes:0.240985):0.240985,Gasterosteus_aculeatus:0.316413):0.05915,Oreochromis_niloticus:0.45):0.08141,(Takifugu_rubripes:0.203847,Tetraodon_nigroviridis:0.224159):1.040705):0.08141,Gadus_morhua:0.16282):0.16282,(Astyanax_mexicanus:0.365376,Danio_rerio:0.365376):0.365376):0.2714825,Lepisosteus_oculatus:0.2714825):0.2714825):0.395016,Petromyzon_marinus:0.790032):0.263344,(Ciona_savignyi:0.8,Ciona_intestinalis:0.8)Cionidae:0.6)Chordata:0.2,(Caenorhabditis_elegans:0.8,Drosophila_melanogaster:0.8):0.8)Coelomata:0.4,Saccharomyces_cerevisiae:1.9)Fungi_Metazoa_group;'

indir = "/u/lafonman/Projects/OrthologsAndParalogs/ConstraintChecker/"
cscheckerbin = "/u/lafonman/Projects/OrthologsAndParalogs/ConstraintChecker/constraint_checker.py"

outdir = "/u/lafonman/data/constraint_checker_results/"


suffixes = ('ultraloose', 'loose', 'default', 'strict', 'ultrastrict')

for treeid in genetreeids:
    
    outfilename = outdir + treeid + ".satcons"
    #------------------
    #this part prevents overriding existing files
    if os.path.isfile(outfilename):
      print outfilename + " exists, skipping"
      continue
    
    
    allexist = True
    for i in range(len(suffixes)):
      if not os.path.isfile(indir + treeid + "-" + suffixes[i] + ".proteinortho-graph"):
	allexist = False
	break
	
    
    if allexist:
      
      
      outfile = open(outfilename, 'w')
      print "WRITING " + outfilename
      print >> outfile, "TREEID=" + treeid
      for i in range(len(suffixes)):
	cmd = "python " + cscheckerbin + ' "' + indir + treeid + "-" + suffixes[i] + '.proteinortho-graph" --speciestree="' + speciesTreeStr + '"'
	
	output = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
	
	outlines = output.split('\n')
	for outline in outlines:
	  if outline.startswith("ISSAT="):
	      print >> outfile, outline.replace("ISSAT", "ISSAT-" + suffixes[i])
	  elif outline.startswith("ISCONS="):
	      print >> outfile, outline.replace("ISCONS", "ISCONS-" + suffixes[i])
	
      
      print "SINGLE CONSTRAINTS DONE"
      
      for i in range(len(suffixes)):
	for j in range(i + 1, len(suffixes)):
	  f1name = indir + treeid + "-" + suffixes[i] + ".proteinortho-graph"
	  f2name = indir + treeid + "-" + suffixes[j] + ".proteinortho-graph"
	  
	  cmd = "python " + cscheckerbin + ' "' + f1name + '" "' + f2name + '" --speciestree="' + speciesTreeStr + '"'
	  output = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
	  print >> outfile, "<" + suffixes[i] + ":" + suffixes[j] + ">"
	  print >> outfile, output
	  print >> outfile, "</" + suffixes[i] + ":" + suffixes[j] + ">"
      outfile.close()
      print "COMBINED CONSTRAINS DONE"
      