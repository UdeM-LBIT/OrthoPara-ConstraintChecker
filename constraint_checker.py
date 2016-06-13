#!/usr/bin/python

import string
import sys
import getopt
import os
from os import listdir
from os.path import isfile, join
import random

from ete2 import TreeNode

"""

Authors: Manuel Lafond, Mark Jones

This script tests satisfiability/consistency with a species tree of a set of orthologies/paralogies/undecided.
The argument --mode={sat,cons,both} specifies which to test.  
If cons or both is set and a species tree is specified, the script will test consistency with the specified species tree.
If cons or both is set and no species tree is satisfied, the script will test whether there exists a species tree with which the orthologies/paralogies are consistent, and will output such a species tree if it exists.
If sat/cons is met, user can output a gene tree meeting the requirements (in ascii and/or newick).

--genes=g1;;g2;;... : list of distinct gene names, all separated by ;;
--speciestree=[newick] : newick string (terminated by ;) having each leaf labeled by a distinct species name
--genespecies=g1:s1;;g2:s2;;... : gene to species mapping.  Each mapping is separated by ';;'.  A mapping has the form GENE_NAME:SPECIES_NAME.  GENE_NAME must be in genes, SPECIES_NAME in a leaf of speciesTree (if speciesTree is specified)
--orthologs=g1:g2;;g1:g3;;g2:g4;;... : all ortholog gene pairs, separated by ':'.  A pair has the form GENE_NAME1:GENE_NAME2
--paralogs=g1:g2;;g1:g3;;g2:g4;;... : all paralog gene pairs, separated by ':'.  A pair has the form GENE_NAME1:GENE_NAME2
--outputascii=[0,1] : will output ascii tree(s) iff set to 1 (default 0)
--outputnewick=[0,1] : will output solution newick for sat and cons iff set to 1 and satisfiability/consistency is met (default 1)
--mode=[sat,cons,both] : If sat, tests satisfiability.  If cons, tests consistency with a species tree (species tree can be specified or unspecified).  If both, tests both. (default both)



Example calls : 
> python constraint_checker.py "--genes=a1;;b1;;c1;;d1;;e1" "--genespecies=a1:a;;b1:b;;c1:c;;d1:d;;e1:e" 
                             "--orthologs=a1:c1;;c1:e1;;e1:d1;;d1:b1" "--paralogs=a1:b1;;a1:d1;;a1:e1;;b1:e1" --outputascii=1 "--speciestree=(c,(d, ((b,a),e)));"
                             
ISSAT=1
((a1,((e1,b1),d1)),c1);

       /-a1
      |
    /dup      /-e1
   |  |    /dup
   |   \spec  \-b1
-spec     |
   |       \-d1
   |
    \-c1
ISCONS=1
(c1,(a1,(d1,(e1,b1))));

    /-c1
-spec
   |   /-a1
    \dup
      |    /-d1
       \spec
          |   /-e1
           \dup
              \-b1


> python constraint_checker.py "--genes=a1;;b1;;c1;;d1;;e1" "--genespecies=a1:a;;b1:b;;c1:c;;d1:d;;e1:e" "--orthologs=a1:c1;;c1:e1;;e1:d1;;d1:b1" "--paralogs=a1:b1;;a1:d1;;a1:e1;;b1:e1" "--speciestree=((a, (b, c)), (d, e));"
ISSAT=1
((a1,((e1,b1),d1)),c1);
ISCONS=0

python constraint_checker.py "--genes=a1;;b1;;c1;;d1;;e1
" "--genespecies=a1:a;;b1:b;;c1:c;;d1:d;;e1:e" "--orthologs=a1:c1;;c1:e1;;e1:d1;;d1:b1" "--paralogs=a1:b1;;a1:d1;;a1:e1;;b1:e1" --outputascii=1 --mode=cons
ISCONS=1
((a1,((e1,b1),d1)),c1);
((a,(e,b),d),c);

       /-a1
      |
    /dup      /-e1
   |  |    /dup
   |   \spec  \-b1
-spec     |
   |       \-d1
   |
    \-c1

      /-a
     |
     |   /-e
   /-|--|
  |  |   \-b
--|  |
  |   \-d
  |
   \-c




This script also has 2 "hidden" modes to take proteinortho graph files as input.  
*** These modes are undocumented ***
*** These modes are intended for calls from a script, such as results_maker.py ***
The graph files contain othology contraints.
If one unnamed argument is passed, the program takes it as the filename of the graph outputted by proteinortho.
In this mode, the script takes all orthologies in the file, and all gene pairs not present are considered paralogies.
No undecided edges in this case.

If 2 unnamed arguments are passed, they are interpreted as 2 graph files.
The program takes all orthologies common to both, all paralogies (ie non-orthologies) common to both, and the rest is left undecided.
The species tree MUST be specified.
A bunch of stats are then outputted.

satisfiability and consistency are tested by default.
TODO: give option to turn this off? Does mode=[sat, cons, both] apply here?

Example calls : 

> python constraint_checker.py ENSGT00390000000020-default.proteinortho-graph
ISSAT=0

> python constraint_checker.py ENSGT00390000000020-default.proteinortho-graph ENSGT00390000000020-loose.proteinortho-graph "--speciestree=((((((((((((((((((((((Pan_troglodytes:0.006667,Homo_sapiens:0.0067):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,Macaca_mulatta:0.037471):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Ictidomys_tridecemlineatus:0.225629):0.01015,Cavia_porcellus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,(((((((Ailuropoda_melanoleuca:0.025614,Mustela_putorius_furo:0.0256):0.0256145,Canis_familiaris:0.051229):0.051229,Felis_catus:0.098612):0.049845,Equus_caballus:0.109397):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508,(((Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,Ovis_aries:0.061796):0.061796):0.025153,Sus_scrofa:0.107275):0.0201675,Vicugna_pacos:0.079):0.0201675):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,((Macropus_eugenii:0.101004,Sarcophilus_harrisii:0.101004):0.021004,Monodelphis_domestica:0.125686):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,(((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,(Ficedula_albicollis:0.085771,Taeniopygia_guttata:0.085771):0.085771):0.199223,Pelodiscus_sinensis:0.489241):0.105143,Anolis_carolinensis:0.4989):0.17):0.149,Xenopus_tropicalis:0.855573):0.155677,Latimeria_chalumnae:0.155677):0.550693,(((((((Xiphophorus_maculatus:0.1204925,Oryzias_latipes:0.240985):0.240985,Gasterosteus_aculeatus:0.316413):0.05915,Oreochromis_niloticus:0.45):0.08141,(Takifugu_rubripes:0.203847,Tetraodon_nigroviridis:0.224159):1.040705):0.08141,Gadus_morhua:0.16282):0.16282,(Astyanax_mexicanus:0.365376,Danio_rerio:0.365376):0.365376):0.2714825,Lepisosteus_oculatus:0.2714825):0.2714825):0.395016,Petromyzon_marinus:0.790032):0.263344,(Ciona_savignyi:0.8,Ciona_intestinalis:0.8)Cionidae:0.6)Chordata:0.2,(Caenorhabditis_elegans:0.8,Drosophila_melanogaster:0.8):0.8)Coelomata:0.4,Saccharomyces_cerevisiae:1.9)Fungi_Metazoa_group;"
...
"""




graphfile1 = ''
graphfile2 = ''
speciesTreeStr = ''
orthologsStr = ''
paralogsStr = ''
genesStr = ''
geneSpeciesStr = ''
outputASCII = False
outputNewick = True
mode = 'both'


#terrible terrible argument parsing, please forgive me
for arg in sys.argv[1:]:
  if arg.startswith("--speciestree="):
    pz = arg.split("=")
    speciesTreeStr = pz[1]
  elif arg.startswith("--genes="):
    pz = arg.split("=")
    genesStr = pz[1]
  elif arg.startswith("--genespecies="):
    pz = arg.split("=")
    geneSpeciesStr = pz[1]
  elif arg.startswith("--mode="):
    pz = arg.split("=")
    mode = pz[1]
  elif arg.startswith("--outputascii="):
    pz = arg.split("=")
    if pz[1] == '1':
      outputASCII = True
  elif arg.startswith("--outputnewick="):
    pz = arg.split("=")
    if pz[1] == '0':
      outputNewick = False
  elif arg.startswith("--orthologs="):
    pz = arg.split("=")
    orthologsStr = pz[1]
  elif arg.startswith("--paralogs="):
    pz = arg.split("=")
    paralogsStr = pz[1]
  elif not arg.startswith('-'):
    if graphfile1 == '':
      graphfile1 = arg
    elif graphfile2 == '':
      graphfile2 = arg

# Map each species name to a leaf.  This avoids using search_nodes, which is slow
# speciesTree is used in every mode, so we build it right here right now
speciesTree = None

if speciesTreeStr != '':
  speciesTree = TreeNode(speciesTreeStr)
  speciesLeavesList = speciesTree.get_leaves()
  speciesLeaves = {}
  for leaf in speciesLeavesList:
    speciesLeaves[leaf.name] = leaf
 


class ConstraintGraph:
  """ Given a set of genes, and and two set of tuples that represent orthology and paralogy relationships between genes,
      builds a graph of orthologies, and a graph of paralogies.  It then becomes possible to build a DS-tree that satisifies
      all these constraints (or detect that it can't be done).
  """
  
  def __init__(self, genes, orthology_relations, paralogy_relations, speciesTree = None, geneSpeciesMapping = None, treelessGeneSpeciesMapping  = None):
    """Constructor.  The 2 graphs are built here, as adjacency lists (self.orthologs and self.paralogs, two dicts with key = gene-string, value = neighbors as a set of gene-strings). 

        :argument genes: The set of genes strings under consideration.  If a gene is in some relation, then it must be in genes.  Example : {"a", "b", "c"}
        :argument orthology_relations: Set of pairs of genes that are orthologous.  Example : {("a", "b"), ("b", "c")}
        :argument paralogy_relations: Same as orthology_relations
        :argument speciesTree: Must be a TreeNode object from ete2.  If one is given, buildDSTree will try to find a DS-tree that is consistent with the speciesTree.
        :argument geneSpeciesMapping: must accompany speciesTree when present.  A dict of gene-string to species tree leaf.
        :argument treelessGeneSpeciesMapping: A dict of gene-string to species name. Used in cases when no speciesTree is given but we want to check consistency. Serves a similar function to geneSpeciesMapping, but genes are mapped to simple species names rather than leaves of a species tree. Used by buildDSTreeAndSpeciesTree.
        """
    self.genes = genes
    self.orthologs = {}                        #key = gene-string, value = set of gene-string neighbors
    self.paralogs = {}                        #same
    self.speciesTree = speciesTree
    self.geneSpeciesMapping = geneSpeciesMapping
    self.treelessGeneSpeciesMapping  = treelessGeneSpeciesMapping 

    
    #build orthology and paralogy neighbors
    for r in orthology_relations:
      if not r[0] in self.orthologs:
        self.orthologs[r[0]] = set()
      if not r[1] in self.orthologs:
        self.orthologs[r[1]] = set()
      self.orthologs[r[0]].add(r[1])
      self.orthologs[r[1]].add(r[0])
    
    
    for r in paralogy_relations:
      if not r[0] in self.paralogs:
        self.paralogs[r[0]] = set()
      if not r[1] in self.paralogs:
        self.paralogs[r[1]] = set()
      self.paralogs[r[0]].add(r[1])
      self.paralogs[r[1]].add(r[0])
    
  def getParalogComponents(self, restrictTo = None):
    """Get connected components of orthology graph, as a list of sets of gene-strings.  
       The idea is that 2 components of the orthology graph are paralogs.
    
        :argument restrictTo: If restrictTo is a set of genes, the components are restricted to this set of genes.
        """
    return self._getComponents(self.orthologs, restrictTo)
    
  def getOrthologComponents(self, restrictTo = None):
    """Get connected components of paralogy graph, as a list of sets of gene-strings.  
    
        :argument restrictTo: If restrictTo is a set of genes, the components are restricted to this set of genes.
        """
    return self._getComponents(self.paralogs, restrictTo)
  
  def _getComponents(self, adjacencies, restrictTo):
    """Get connected components of graph corresponding to the given adjacencies, as a list of sets of gene-strings.  
   
        :argument restrictTo: If restrictTo is a set of genes, the components are restricted to this set of genes.
        """
    components = []
    
    #We shall perform a DFS on the graph.
    if restrictTo is None:
      unvisited = self.genes.copy()
    else:
      unvisited = restrictTo.copy()
    
    while len(unvisited) > 0:
      v = unvisited.pop()
      curcomponent = set()
      curcomponent.add(v)
      self._doDFS(adjacencies, v, unvisited, curcomponent, restrictTo)
      components.append(curcomponent)
      
    return components
      
  def _doDFS(self, adjacencies, v, unvisited, curcomponent, restrictTo):
    """Recursive function for a DFS.  
    
        :argument adjacencies: The graph
        :argument v: The node we're currently visiting.
        :argument unvisited: Nodes yet to visit in the DFS
        :argument curComponent: The component we're currently filling
        :restrictTo: see above
        """
    if not v in adjacencies:
      return
    
    for v2 in adjacencies[v]:
      if restrictTo == None or v2 in restrictTo:
        if v2 in unvisited:
          curcomponent.add(v2)
          unvisited.remove(v2)
          self._doDFS(adjacencies, v2, unvisited, curcomponent, restrictTo)

  
  def buildDSTree(self):
    """Build one valid DS-tree for the relationships given (and possibly the speciesTree) and return it.
      If none can be built, returns False instead.
      The DS-tree built prioritizes dups first.
        """
    dstree = TreeNode()
    
    startcc = self.genes.copy()
    hasPassed = self._buildDSTree(startcc, dstree)
    
    #in case you are wondering, what's below is (hasPassed ? dstree : False)
    return (dstree if hasPassed else False)


  def buildDSAndSpeciesTree(self):
    """Build one valid DS-tree for the relationships given, and a species tree with which it is consistent, and returns both trees.
      If no such pair of trees can be built, returns False instead.
      The DS-tree built prioritizes dups first.
        """
    
    startcc = self.genes.copy()
    dstree = TreeNode()
    dstree.set = set(startcc)
    dsNodeList =[dstree]
    constructedSpeciesTree = TreeNode()	
    constructedSpeciesTree.set = set(self.treelessGeneSpeciesMapping[x] for x in startcc)
    hasPassed = self._buildDSAndSpeciesTree(startcc, dsNodeList, constructedSpeciesTree)

    #Return dstree and species tree as an ordered pair
    treepair = [dstree, constructedSpeciesTree]
    
    #in case you are wondering, what's below is (hasPassed ? treepair : False)
    return (treepair if hasPassed else False)

  
  def _buildDSTree(self, curCC, curNode):
    """Recursively build the DS-tree
    
        :argument curCC: The subgraph we're working on.  It's a list of genes.
        :argument curNode: The node in the DS-tree we're working under.  We'll add children to this node here.
        """
        
        
    if len(curCC) == 1:
      c = curCC.pop()
      curNode.name = c
      return True
    
    curtype = "dup"
    newComponents = self.getParalogComponents(curCC)
    if len(newComponents) == 1:
      curtype = "spec"
      newComponents = self.getOrthologComponents(curCC)
      
      #If the speciesTree is given, we might not want to split the components
      #by speciation maximally.  Some components might need to be merged together
      #as they can't be joined by speciation, and we hope we'll be able to join them by duplication later.
      #_getSpeciationPartition takes care of this merging.
      if self.speciesTree != None:
        newComponents = self._getSpeciationPartition(newComponents)
      
      #Orthology + Paralogy graphs are connected and non-trivial => no solution
      if len(newComponents) == 1:
        return False
    
    curNode.name = curtype
    hasPassed = True
    #TODO : we might want to stop as soon as hasPassed is set to false
    for cc in newComponents:
      child = curNode.add_child()
      hasPassed = self._buildDSTree(cc, child) and hasPassed
      
    return hasPassed
      

  def _buildDSAndSpeciesTree(self, curCC, curDSNodes, curSpeciesNode):
    """Recursively build the DS-tree and species tree
    
        :argument curCC: The subgraph we're working on.  It's a list of genes. TODO cut this?
        :argument curDSNodes: The set of nodes in the DS-tree we're working under.  We'll add children to these nodes here.
        :argument curSpeciesNode: The node in the species tree we're working under. We'll add children to this node here.
        """

    # Unlike _BuildSpeciesTree, where we have a single dstree node, here we have a set of nodes to work with.

    # Firstly, detect duplication events that happen before the next split in the species tree.
    # This is done by finding nodes in curDSNodes whose gene sets can be separated without breaking any orthology relations.
    # For such a node, we turn iot into a duplication node and split its genes up among its children.
    newDSNodeList = [] 
    for node in curDSNodes:   
      #TODO hopefully will not need to cast node.set as a set, once construction of node.set elsewhere is fixed up.
      newComponents = self.getParalogComponents(set(node.set))
      if len(newComponents) > 1: # node is a duplication node; find its children, remove it from node list and add its children to list.
        node.name = 'dup'
        for cc in newComponents:
          child = node.add_child()
          child.set = set(cc)
          newDSNodeList.append(child)
      else:					#node is not a leaf or duplication node; keep it in the node list.
         newDSNodeList.append(node)
    curDSNodes = newDSNodeList

    # Now that the duplication nodes are out of the way, it's time to find a speciation partition!

    # If we only have one species left, then as genes from the same species must be paralogs, we must have that all genes are isolated.
    #Note this relies on the fact that we already reduced the nodes into their smallest possible connected components, earlier in this method.
    if len(curSpeciesNode.set) == 1:
      hasPassed = True
      for node in curDSNodes:
        if len(node.set) != 1:
          hasPassed = False
      # If there is no problem, complete the DStree and species tree at these nodes
      if hasPassed:
        curSpeciesNode.name = set(curSpeciesNode.set).pop()
        for node in curDSNodes:
          #Use set(node.set).pop() rather than node.set.pop() here because pop() deletes the item from the set, which can cause wierd errors. I think I tracked down the bug anyway, but still removing the item from the set isn't really what we're meant to be doing.
          node.name = set(node.set).pop()  
          #node.name = node.set.pop()
        return True
      else:
        return False

    # If we have multiple species, we now find a speciation partition and recurse
    else:
      #Find the 'disconnected components' (i.e. components in the complement) for each node set.
      disconnectedComponents = []
      for node in curDSNodes:
        tempComponents = self.getOrthologComponents(node.set)
        disconnectedComponents.extend(tempComponents)
  
      # find a speciation partition based on disconnected components
      speciationComponents = self._getSpeciationPartitionForUnknownTree(disconnectedComponents)
      # If the speciation partition has 1 component, there is no valid way of splitting up species, and we have no solution.
      if len(speciationComponents) == 1:
        return False

     #identify speciation nodes based on this partition, and recurse on each set of species
      hasPassed = True
      for sc in speciationComponents:
        newDSList = []
        speciesChild = curSpeciesNode.add_child()
        speciesChild.set = set(self.treelessGeneSpeciesMapping[x] for x in sc)
        for node in curDSNodes:
          intersection = node.set.intersection(sc)
          #NOTE code is a bit fragile here: need len(node.set) > 0 because otherwise a node can get passed down through the recursion and come back up empty, and will then get sent places it shouldn't.
          if intersection == node.set and len(node.set) > 0: #do nothing with this node for now but keep it in the node list.
            newDSList.append(node)
          else:
            if len(intersection) > 0: #node.set has some (but not all) genes within speciation component; process these genes for this speciation component
              node.name = 'spec'
              child = node.add_child()
              child.set = intersection
              newDSList.append(child)
        #recurse on a new instance, restricted to the genes in sc
        hasPassed = self._buildDSAndSpeciesTree(sc, newDSList, speciesChild) and hasPassed
      return hasPassed




  def _getSpeciationPartition(self, newComponents):
    """Partitions the components into "unrelated" components - that is, 
       components X and Y are 'merged' iff lca(s(X)) is related to lca(s(Y))
    
        :argument newComponents: The components we're working with.  Name comes from historical facts and lazyness to change it.  TODO !
        """
    newerComponents = []        #the partition
    
    #Here we build a dict with key = speciesTree node, value = list of components X s.t. lca(s(X)) = that node
    speciesComponents = {}
    for c in newComponents:
      s = self.getSpeciesLCA(c)
      if s not in speciesComponents:
        speciesComponents[s] = []
      speciesComponents[s].append(c)
    
    #a BFS will find the partitions
    self._doBFSForHighestSpecies(self.speciesTree, speciesComponents, newerComponents, None)
    
    return newerComponents
      
    
  def getSpeciesLCA(self, genes):
    """Returns lca(s(genes))
    
        :argument genes: The genes set
        """
    s = set()
    for g in genes:
      s.add( self.geneSpeciesMapping[g] )

    # BUG FIX: node.get_common_ancestor(s) (below) has either changed or does not work as intended.
    # Intended outcome: return lca of node and s.
    # Actual outcome: return lca of s.
    # For this reason, am changing code so that we do not remove node from the set s when we get it.
    # 29/04/2016 Mark Jones
    node = set(s).pop()
    #node = s.pop()    

    if len(s) == 0:
      return node
    else:
      return node.get_common_ancestor( s )
    
  def _doBFSForHighestSpecies(self, curNode, speciesComponents, returnedComponents, curComponent = None):
    """Recursive BFS function to find a speciation partition.
       We identify a new component when no node above curNode was in speciesComponents, but curNode is.
       Then, all nodes under curNode that are in speciesComponents will be part of this new component.
       See it as a maximal subtree s.t. the root is in speciesComponents
    
        :argument curNode: Current node in the BFS
        :argument speciesComponents: A dict with key = species tree node, value = list of initial components having lca on that node
        :argument returnedComponents: We fill this thing with all the components we create
        :argument curComponent:component we're working on (the new component identified as in the function description, if any).
        """
    
    #TODO : here we set curComponent as it was when entering before exiting function call.  Is that necessary ???
    curComponentBefore = curComponent
    
    curComponentWasCreatedHere = False
    
    if curNode in speciesComponents:
      
      #hey, we'll start a new component here
      if curComponent is None:
        curComponent = set()
        curComponentWasCreatedHere = True
        
      #merge speciesComponents with current component
      for c in speciesComponents[curNode]:
        curComponent.update(c)
    
    
    for child in curNode.children:
      self._doBFSForHighestSpecies(child, speciesComponents, returnedComponents, curComponent)
      
    if curComponentWasCreatedHere:
      returnedComponents.append(curComponent)
      
    curComponent = curComponentBefore
      
       
  def _getSpeciationPartitionForUnknownTree(self, initialComponents):
    """A version of _getSpeciationPartition when we don't have a known species tree to work with, and have to construct one for ourselves.
  #Given set of components over a set of genes, return a partition of the same set of genes, such that every component is contained within a single part, and no two parts have genes from the same species.
      """
    newerComponents = []	#the partition, initially set to be same as set of initial components
    for x in initialComponents:
      newerComponents.append(x) 

    #now merge components if they have genes from the same species
    # possibly terrible code, I am sorry!
    for dc in initialComponents:
      newPart = []
      deletionSet = []
      for dc2 in newerComponents:
        merge = False
        for x in dc:
          for y in dc2:
            if self.treelessGeneSpeciesMapping[x] == self.treelessGeneSpeciesMapping[y]:
              merge = True
              break
          if merge == True:
            break
        if merge == True:
          newPart.extend(dc2)
          deletionSet.append(dc2)
       #update newerComponents, with the 'deleted' parts merged into newPart
      if len(newPart) >0:
        for dc2 in deletionSet:
          newerComponents.remove(dc2)
        newerComponents.append(newPart)
    return newerComponents


def getGeneSpeciesMapping(genes, speciesTree, separator = '__', speciesLeaves = {}):
  """ Returns a map in which the keys are the gene names, and the value its corresponding species.
      Values in genes are expected to have the format [GENENAME][SEPARATOR][SPECIESTREE]
    
        :argument genes: List of gene names that each include separator and species
        :argument speciesTree: a TreeClass object
        :argument separator: [GENENAME][SEPARATOR][SPECIESTREE]
        """
  
  mapping = {}

  for g in genes:
    if separator is None:
      genename = g
      speciesname = g
    else:
      pz = g.split(separator)
      genename = pz[0]
      speciesname = pz[1]
    #search_nodes is slooooow, avoid it whenever possible

    if speciesname in speciesLeaves:
      mapping[genename] = speciesLeaves[speciesname]
    else:
      mapping[genename] = speciesTree.search_nodes(name=speciesname)[0]  
  
  return mapping
    





def getTreelessGeneSpeciesMapping(genes, separator = '__', speciesLeaves = {}):
  """ Returns a map in which the keys are the gene names, and the value its corresponding species.
      Values in genes are expected to have the format [GENENAME][SEPARATOR][SPECIESNAME]
      Unlike getGeneSpeciesMapping, in which we are given a specie tree and genes are mapped to leaves of the tree, here we do not have a species tree and genes are simply mapped to species names.
    
        :argument genes: List of gene names that each include separator and species
        :argument separator: [GENENAME][SEPARATOR][SPECIESTREE]
        """
  
  mapping = {}

  for g in genes:
    if separator is None:
      genename = g
      speciesname = g
    else:
      pz = g.split(separator)
      genename = pz[0]
      speciesname = pz[1]
    mapping[genename] = speciesname
  return mapping


#TODO make this method generate instances without species trees?
def generateRandomProblem(nbGenes, nbSpecies, orthologProb = 0.5, paralogProb = 0.4):
  """ Generates a random set of genes, orthologs, paralogs and species tree, ready to be input in ConstraintGraph class.
      Returned value has the form 
      { "genes" : geneset, "orthologs" : orthologs, "paralogs" : paralogs, "speciesTree" : speciesTree }
      geneset items have the form [GENENAME]:[SPECIESNAME]
      Gene species are attributed randomly, though each species has at least one gene.
      If two genes have the same species, they end up in paralogs, always.
      
        :argument nbGenes: Number of genes to generate
        :argument orthologProb: chances for 2 genes to be a pair in orthologs
        :argument paralogProb: chances for 2 genes to be a pair in paralogs
        """
  speciesnames = range(0, nbSpecies)
  speciesTree = TreeNode()
  speciesTree.populate(nbSpecies, speciesnames)
  
  for node in speciesTree:
    if not node.name is None:
      node.name = str(node.name)
  
  genes = []
  
  #first add one gene per species
  for i in range(nbSpecies):
    genes.append("g" + str(len(genes)) + ":" + str(i))
    
  #then fill in the rest with random species genes
  paralogs = set()
  orthologs = set()
  
  for i in range(nbGenes - nbSpecies):
    s = random.randint(0, nbSpecies - 1)
    genes.append("g" + str(len(genes)) + ":" + str(s))
    
  #and here we decide of random relationships
  paralogProb += orthologProb 
  for i in range(nbGenes):
    g1 = genes[i]
    pz = g1.split(":")
    g1name = pz[0]
    g1species = pz[1]
    
    for j in range(i + 1, nbGenes):
      g2 = genes[j]
      px = g2.split(":")
      g2name = px[0]
      g2species = px[1]
      
      if g1species == g2species:
        paralogs.add( (g1, g2) )
      else:
        p = random.random()
        
        if p < orthologProb:
          orthologs.add( (g1, g2) )
        elif p >= orthologProb and p < paralogProb:
          paralogs.add( (g1, g2) )

  geneset = set()
  geneset.update(genes)
  return { "genes" : geneset, "orthologs" : orthologs, "paralogs" : paralogs, "speciesTree" : speciesTree }

def loadGraphFile(graphfile):
    """ Loads a graph file outputted by proteinortho.  Mainly does line parsing to extract ortholog gene pairs.
        Return value is a map of the form 
        { "genes" : genes, "orthologs" : orthologs, "paralogs" : paralogs, "geneSpeciesMapping" : geneSpeciesMapping, "speciesGenes" : speciesGenes }
        genes is a set of genes (not including species name)
        geneSpeciesMapping maps gene names to species names
        speciesGenes maps species to the set of genes they contain
        :argument graphfile: filename
        """
    genes = set()
  
  
    
    gfile = open(graphfile)
    
    orthologs = set()
    paralogs = set()
    
    curspecies1 = ''
    curspecies2 = ''
    
    geneSpeciesMapping = {}
    speciesGenes = {}
    
    for line in gfile:
      line = line.replace('\n', '')
        
      
      if line.startswith('#'):
        line = line.replace('#', '').strip().replace('.fasta', '')
        pz = line.split('\t')
        curspecies1 = pz[0]
        curspecies2 = pz[1]
        
        if curspecies1 not in speciesGenes:
          speciesGenes[curspecies1] = set()
        if curspecies2 not in speciesGenes:
          speciesGenes[curspecies2] = set()
          
      elif curspecies1 != '' and curspecies2 != '':
        pz = line.split('\t')
        if pz[0] not in genes:
          genes.add(pz[0])
        if pz[1] not in genes:
          genes.add(pz[1])
          
        if not (pz[0], pz[1]) in orthologs and not (pz[1], pz[0]) in orthologs:
          orthologs.add( (pz[0], pz[1]) )
        geneSpeciesMapping[pz[0]] = curspecies1
        geneSpeciesMapping[pz[1]] = curspecies2
        
        if not pz[0] in speciesGenes[curspecies1]:
          speciesGenes[curspecies1].add(pz[0])
        if not pz[1] in speciesGenes[curspecies2]:
          speciesGenes[curspecies2].add(pz[1])
        
        
    #add forced paralogs (those of the same species
    for s in speciesGenes:
      
      for g1 in speciesGenes[s]:
        for g2 in speciesGenes[s]:
          if g1 != g2:
            if not (g1, g2) in paralogs and not (g2, g1) in paralogs:
              paralogs.add( (g1, g2) )
            
    #then, what's not ortholog is considered paralog
    for g1 in genes:
      for g2 in genes:
        if g1 != g2:
          if (g1, g2) not in orthologs and (g2, g1) not in orthologs and not (g1, g2) in paralogs and not (g2, g1) in paralogs:
            paralogs.add( (g1, g2) )
            

    gfile.close()
    
    return { "genes" : genes, "orthologs" : orthologs, "paralogs" : paralogs, "geneSpeciesMapping" : geneSpeciesMapping, "speciesGenes" : speciesGenes }


#--------------------------------------------------------------------------------------------------------
# STUDY A SINGLE GRAPH AND OUTPUT WHETHER IT IS SATISFIABLE, AND CONSISTENT WITH THE SPECIES TREE
#--------------------------------------------------------------------------------------------------------
if graphfile1 != '' and graphfile1 != 'random' and graphfile2 == '':
  
  gp = loadGraphFile(graphfile1)
  
  graph = ConstraintGraph(gp['genes'], gp['orthologs'], gp['paralogs'])
  t = graph.buildDSTree()
  
  #TODO : we just output whether tree is ok or not
  if t == False:
    print "ISSAT=0"
  else:
    print "ISSAT=1"
    
  #TODO : THIS CODE IS REPEATED IN THE CASE OF 2 GRAPH FILES
    
    #in case we can't map a gene, we will remove it
  g1mapping = gp['geneSpeciesMapping']
  geneNamesToErase = set()
  geneSpeciesMapping = {}

    # Only construct geneSpeciesMapping if species tree exists
  if speciesTree != None:
    for g in gp['genes']:
      if g1mapping[g] not in speciesLeaves:
        geneNamesToErase.add(g)
      else:
        geneSpeciesMapping[g] = speciesLeaves[g1mapping[g]]

    genes = gp['genes']
    for gerase in geneNamesToErase:
      genes.remove(gerase)
  
    graph = ConstraintGraph(gp['genes'], gp['orthologs'], gp['paralogs'], speciesTree, geneSpeciesMapping)
    t = graph.buildDSTree()
    if t == False:
      print "ISCONS=0"
    else:
      print "ISCONS=1"
      print "DSTREECONS=" + t.write()   
        
    #If no species tree exists, we simply use the mapping given by g1mapping as our "treeless" gene species mapping
  if speciesTree == None:
    #If no species tree exists, we simply use the mapping given by g1mapping as our "treeless" gene species mapping
    treelessGeneSpeciesMapping = g1mapping
    graph = ConstraintGraph(gp['genes'], gp['orthologs'], gp['paralogs'], None, None, treelessGeneSpeciesMapping)
    treepair = graph.buildDSAndSpeciesTree()
    if treepair != False:
      t = treepair[0]
      s = treepair[1]
      print "ISCONS=1"
      print "DSTREECONS=" + t.write()
      print "SPECTREECONS=" + s.write()     

	  #below = code to check whether relations also agree with some fixed, hard-coded species tree
      #speciesTreeStrTest = '((((((((((((((((((((((Pan_troglodytes:0.006667,Homo_sapiens:0.0067):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,Macaca_mulatta:0.037471):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Ictidomys_tridecemlineatus:0.225629):0.01015,Cavia_porcellus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,(((((((Ailuropoda_melanoleuca:0.025614,Mustela_putorius_furo:0.0256):0.0256145,Canis_familiaris:0.051229):0.051229,Felis_catus:0.098612):0.049845,Equus_caballus:0.109397):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508,(((Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,Ovis_aries:0.061796):0.061796):0.025153,Sus_scrofa:0.107275):0.0201675,Vicugna_pacos:0.079):0.0201675):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,((Macropus_eugenii:0.101004,Sarcophilus_harrisii:0.101004):0.021004,Monodelphis_domestica:0.125686):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,(((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,(Ficedula_albicollis:0.085771,Taeniopygia_guttata:0.085771):0.085771):0.199223,Pelodiscus_sinensis:0.489241):0.105143,Anolis_carolinensis:0.4989):0.17):0.149,Xenopus_tropicalis:0.855573):0.155677,Latimeria_chalumnae:0.155677):0.550693,(((((((Xiphophorus_maculatus:0.1204925,Oryzias_latipes:0.240985):0.240985,Gasterosteus_aculeatus:0.316413):0.05915,Oreochromis_niloticus:0.45):0.08141,(Takifugu_rubripes:0.203847,Tetraodon_nigroviridis:0.224159):1.040705):0.08141,Gadus_morhua:0.16282):0.16282,(Astyanax_mexicanus:0.365376,Danio_rerio:0.365376):0.365376):0.2714825,Lepisosteus_oculatus:0.2714825):0.2714825):0.395016,Petromyzon_marinus:0.790032):0.263344,(Ciona_savignyi:0.8,Ciona_intestinalis:0.8):0.6):0.2,(Caenorhabditis_elegans:0.8,Drosophila_melanogaster:0.8):0.8):0.4,Saccharomyces_cerevisiae:1.9);'
      #speciesTreeTest = TreeNode(speciesTreeStrTest)
      #speciesLeavesListTest = speciesTreeTest.get_leaves()
      #speciesLeavesTest = {}
      #for leaf in speciesLeavesListTest:
      #  speciesLeavesTest[leaf.name] = leaf
      #geneNamesToEraseTest = set()
      #geneSpeciesMappingTest = {}
      #for g in gp['genes']:
      #  if g1mapping[g] not in speciesLeavesTest:
      #    geneNamesToEraseTest.add(g)
      #  else:
      #    geneSpeciesMappingTest[g] = speciesLeavesTest[g1mapping[g]]

      #genes = gp['genes']
      
      #for gerase in geneNamesToEraseTest:
      #  genes.remove(gerase)
                
      #graphTest = ConstraintGraph(gp['genes'], gp['orthologs'], gp['paralogs'], speciesTreeTest, geneSpeciesMappingTest)
      #q = graphTest.buildDSTree()
      #if q == False:
      #  print "SPECTREECONSNOST=" + s.write()
      
        
    else:
      
      print "ISCONS=0"

   
#--------------------------------------------------------------------------------------------------------
# WE'VE GOT TWO GRAPH FILES, AND WE JUST KEEP ORTHOLOGS/PARALOGS THEY HAVE IN COMMON
# THEN WE OUTPUT A BUNCH OF STATS
#--------------------------------------------------------------------------------------------------------
elif graphfile1 != '' and graphfile2 != '':
  gp1 = loadGraphFile(graphfile1)
  gp2 = loadGraphFile(graphfile2)
  
  genes = gp1['genes'].copy()
  orthologs = set()
  paralogs = set()
  
  #orthologs are those found in both graphs
  strortho = ''
  strpara = ''
  for ort1 in gp1['orthologs']:
    if ((ort1[0], ort1[1]) in gp2['orthologs'] or (ort1[1], ort1[0]) in gp2['orthologs']) and (ort1[0], ort1[1]) not in orthologs and (ort1[1], ort1[0]) not in orthologs:
      orthologs.add( ort1 )
      if strortho != '':
        strortho += ','
      strortho += ort1[0] + ':' + ort1[1]
  
  #same with paralogs
  for par1 in gp1['paralogs']:
    if ((par1[0], par1[1]) in gp2['paralogs'] or (par1[1], par1[0]) in gp2['paralogs']) and (par1[0], par1[1]) not in paralogs and (par1[1], par1[0]) not in paralogs:
      paralogs.add( par1 )
      if strpara != '':
        strpara += ','
      strpara += par1[0] + ':' + par1[1]


  #in case we can't map a gene, we will remove it
  g1mapping = gp1['geneSpeciesMapping']
  geneNamesToErase = set()
  
  for g in genes:
    if g not in gp2['genes']:
      geneNamesToErase.add(g)
  
  # Only construct geneSpeciesMapping if species tree exists
  geneSpeciesMapping = {}
  if speciesTree != None:
    for g in genes:
      if g1mapping[g] not in speciesLeaves:
        geneNamesToErase.add(g)
      else:
        geneSpeciesMapping[g] = speciesLeaves[g1mapping[g]]


  # If no species tree exists, we simply use the mapping given by g1mapping as our "treeless" gene species mapping
  if speciesTree == None:
    geneSpeciesStrItems = geneSpeciesStr.split(';;')
    treelessGeneSpeciesMapping = g1mapping

  
  orthologsToErase = set()
  paralogsToErase = set()

  for gerase in geneNamesToErase:
    genes.remove(gerase)
    for ort in orthologs:
      if ort[0] == gerase or ort[1] == gerase:
        orthologsToErase.add(ort)
    for par in paralogs:
      if par[0] == gerase or par[1] == gerase:
        paralogsToErase.add(par)

  
  for oerase in orthologsToErase:
    orthologs.remove(oerase)
  for perase in paralogsToErase:
    paralogs.remove(perase)

  

  print "NBGENES=" + str(len(genes))
  print "NBGENES1=" + str(len(gp1['genes']))
  print "NBGENES2=" + str(len(gp2['genes']))
  print "NBSPECIES1=" + str(len(gp1['speciesGenes']))
  print "NBSPECIES2=" + str(len(gp2['speciesGenes']))
  print "NBORTHOLOGS1=" + str(len(gp1['orthologs']))
  print "NBORTHOLOGS2=" + str(len(gp2['orthologs']))
  print "NBPARALOGS1=" + str(len(gp1['paralogs']))
  print "NBPARALOGS2=" + str(len(gp2['paralogs']))
  print "NBORTHOLOGS=" + str(len(orthologs))
  print "NBPARALOGS=" + str(len(paralogs))
  print "ORTHOLOGS=" + strortho
  print "PARALOGS=" + strpara
  
  
  
  graph = ConstraintGraph(genes, orthologs, paralogs)
  t = graph.buildDSTree()
    
  if t == False:
    print "ISSAT=0"
  else:
    print "ISSAT=1"
    print "DSTREE=" + t.write()

  
  if speciesTree == None:
      graph = ConstraintGraph(genes, orthologs, paralogs, None, None, treelessGeneSpeciesMapping)
      treepair = graph.buildDSAndSpeciesTree()
      if treepair != False:
        t = treepair[0]
        s = treepair[1]
        print "ISCONS=1"
        print "DSTREECONS=" + t.write()
        print "SPECTREECONS=" + s.write()
        
        #speciesTreeStrTest = '((((((((((((((((((((((Pan_troglodytes:0.006667,Homo_sapiens:0.0067):0.00225,Gorilla_gorilla:0.008825):0.00968,Pongo_abelii:0.018318):0.00717,Nomascus_leucogenys:0.025488):0.00717,Macaca_mulatta:0.037471):0.021965,Callithrix_jacchus:0.066131):0.05759,Tarsius_syrichta:0.137823):0.011062,(Microcebus_murinus:0.092749,Otolemur_garnettii:0.129725):0.035463):0.015494,Tupaia_belangeri:0.186203):0.004937,(((((Mus_musculus:0.084509,Rattus_norvegicus:0.091589):0.197773,Dipodomys_ordii:0.211609):0.022992,Ictidomys_tridecemlineatus:0.225629):0.01015,Cavia_porcellus:0.148468):0.025746,(Oryctolagus_cuniculus:0.114227,Ochotona_princeps:0.201069):0.101463):0.015313):0.020593,(((((((Ailuropoda_melanoleuca:0.025614,Mustela_putorius_furo:0.0256):0.0256145,Canis_familiaris:0.051229):0.051229,Felis_catus:0.098612):0.049845,Equus_caballus:0.109397):0.006219,(Myotis_lucifugus:0.14254,Pteropus_vampyrus:0.113399):0.033706):0.004508,(((Tursiops_truncatus:0.064688,(Bos_taurus:0.061796,Ovis_aries:0.061796):0.061796):0.025153,Sus_scrofa:0.107275):0.0201675,Vicugna_pacos:0.079):0.0201675):0.011671,(Erinaceus_europaeus:0.221785,Sorex_araneus:0.269562):0.056393):0.021227):0.023664,(((Loxodonta_africana:0.082242,Procavia_capensis:0.155358):0.02699,Echinops_telfairi:0.245936):0.049697,(Dasypus_novemcinctus:0.116664,Choloepus_hoffmanni:0.096357):0.053145):0.006717):0.234728,((Macropus_eugenii:0.101004,Sarcophilus_harrisii:0.101004):0.021004,Monodelphis_domestica:0.125686):0.2151):0.071664,Ornithorhynchus_anatinus:0.456592):0.109504,(((((Gallus_gallus:0.041384,Meleagris_gallopavo:0.041384):0.041384,Anas_platyrhynchos:0.082768):0.082768,(Ficedula_albicollis:0.085771,Taeniopygia_guttata:0.085771):0.085771):0.199223,Pelodiscus_sinensis:0.489241):0.105143,Anolis_carolinensis:0.4989):0.17):0.149,Xenopus_tropicalis:0.855573):0.155677,Latimeria_chalumnae:0.155677):0.550693,(((((((Xiphophorus_maculatus:0.1204925,Oryzias_latipes:0.240985):0.240985,Gasterosteus_aculeatus:0.316413):0.05915,Oreochromis_niloticus:0.45):0.08141,(Takifugu_rubripes:0.203847,Tetraodon_nigroviridis:0.224159):1.040705):0.08141,Gadus_morhua:0.16282):0.16282,(Astyanax_mexicanus:0.365376,Danio_rerio:0.365376):0.365376):0.2714825,Lepisosteus_oculatus:0.2714825):0.2714825):0.395016,Petromyzon_marinus:0.790032):0.263344,(Ciona_savignyi:0.8,Ciona_intestinalis:0.8):0.6):0.2,(Caenorhabditis_elegans:0.8,Drosophila_melanogaster:0.8):0.8):0.4,Saccharomyces_cerevisiae:1.9);'
        #speciesTreeTest = TreeNode(speciesTreeStrTest)

        #speciesLeavesListTest = speciesTreeTest.get_leaves()
        #speciesLeavesTest = {}        
        #for leaf in speciesLeavesListTest:
        #  speciesLeavesTest[leaf.name] = leaf
        #geneNamesToEraseTest = set()
        #geneSpeciesMappingTest = {}
        #for g in genes:
        #  if g1mapping[g] not in speciesLeavesTest:
        #    geneNamesToEraseTest.add(g)
        #  else:
        #    geneSpeciesMappingTest[g] = speciesLeavesTest[g1mapping[g]]

        #orthologsToEraseTest = set()
        #paralogsToEraseTest = set()
        

        #for gerase in geneNamesToEraseTest:
        #  genes.remove(gerase)
        #  for ort in orthologs:
        #    if ort[0] == gerase or ort[1] == gerase:
        #      orthologsToEraseTest.add(ort)
        #  for par in paralogs:
        #    if par[0] == gerase or par[1] == gerase:
        #      paralogsToEraseTest.add(par)

  
        #for oerase in orthologsToEraseTest:
        #  orthologs.remove(oerase)
        #for perase in paralogsToEraseTest:
        #  paralogs.remove(perase)    
          
        #graphTest = ConstraintGraph(genes, orthologs, paralogs, speciesTreeTest, geneSpeciesMappingTest)

        #q = graphTest.buildDSTree()
        #if q == False:
        #  print "SPECTREECONSNOST=" + s.write()   
        #if q == False:
        #  print "SPECTREECONSNOSTERROR" 
        #else:
        #  print "SPECTREECONSNOSTOK" 

      else:
        print "ISCONS=0"

  if speciesTree != None:    
    
    graph = ConstraintGraph(genes, orthologs, paralogs, speciesTree, geneSpeciesMapping)
    t = graph.buildDSTree()
    
    if t == False:
      print "ISCONS=0"
    else:
      print "ISCONS=1"
      print "DSTREECONS=" + t.write()

if graphfile1 == '' and graphfile2 == '':
  if genesStr == '':
    print "Argument error.  No genes specified."
    sys.exit()
  
  
  genesList = genesStr.split(';;')
  genes = set()
  genes.update(genesList)
  
  #PARSE SPECIES TREE AND MAPPING
  
  if speciesTree != None and geneSpeciesStr != '':
    #speciesTree = TreeNode(speciesTreeStr)   already done way above here
    geneSpeciesStrItems = geneSpeciesStr.split(';;')
   
    geneSpeciesMapping = getGeneSpeciesMapping(geneSpeciesStrItems, speciesTree, ":")

  # Create gene mapping for the cases when we are given a geneSpeciesStr but no speciesTree (i.e. the cases where we want to reconstruct the species tree)
  if speciesTree == None and geneSpeciesStr != '':
    geneSpeciesStrItems = geneSpeciesStr.split(';;')
    treelessGeneSpeciesMapping = getTreelessGeneSpeciesMapping(geneSpeciesStrItems, ":")
    
  #PARSE ORTHOLOGS AND PARALOGS
  orttmp = orthologsStr.split(';;')
  partmp = paralogsStr.split(';;')
  
  orthologs = set()
  paralogs = set()
  
  #TODO : error verification
  for otmp in orttmp:
    pz = otmp.split(':')
    orthologs.add( (pz[0], pz[1] ) )
  for ptmp in partmp:
    pz = ptmp.split(':')
    paralogs.add( (pz[0], pz[1] ) )
    
  #DO SATISFIABILITY  
  if mode == 'sat' or mode == 'both':
    graph = ConstraintGraph(genes, orthologs, paralogs)
    t = graph.buildDSTree()

    
    if t != False:
      print "ISSAT=1"
      
      if outputNewick:
        print t.write(format=9)
      
      if outputASCII:
        print t.get_ascii(show_internal = True)
    else:
      print "ISSAT=0"
      
      
  #DO CONSISTENCY
  if (mode == 'cons' or mode == 'both'):
    # if no species tree, try and build one in addition to DS-tree
    if speciesTree == None:
      graph = ConstraintGraph(genes, orthologs, paralogs, None, None, treelessGeneSpeciesMapping)
      treepair = graph.buildDSAndSpeciesTree()
      if treepair != False:
        t = treepair[0]
        s = treepair[1]
        print "ISCONS=1"
      
        if outputNewick:
          print t.write(format=9)
          print s.write(format=9)
      
        if outputASCII:
          print t.get_ascii(show_internal = True)
          print s.get_ascii(show_internal = True)
      else:
        print "ISCONS=0"

    else:
      graph = ConstraintGraph(genes, orthologs, paralogs, speciesTree, geneSpeciesMapping)
      t = graph.buildDSTree()

    
      if t != False:
        print "ISCONS=1"
      
        if outputNewick:
          print t.write(format=9)
      
        if outputASCII:
          print t.get_ascii(show_internal = True)
      else:
        print "ISCONS=0"
        
#--------------------------------------------------------------------------------------------------------
# PERFORM RANDOM TESTS
# THIS PROBABLY DOESN'T WORK ANYMORE
#--------------------------------------------------------------------------------------------------------


'''
if graphfile1 == 'random':

    print "P(Ort),P(Par),nbG,nbS,Satisfiable,Consistent"

    probs = [ 0.01, 0.05, 0.1, 0.2, 0.3, 0.4 ]
    spgenefactor = 0.75
    nbIters = 20

    for p_o in range(len(probs)):
      for p_p in range(len(probs)):
        prob_orth = probs[p_o]
        prob_para = probs[p_p]
        for nbGenes in range(4, 100):
          for it in range(nbIters):
            nbSpecies = int(nbGenes * spgenefactor)
            params = generateRandomProblem(nbGenes, nbSpecies, prob_orth, prob_para)
            
            genes = params["genes"]
            speciesTree = params["speciesTree"]
            orthologs = params["orthologs"]
            paralogs = params["paralogs"]
        
            graph = ConstraintGraph(genes, orthologs, paralogs)
            t = graph.buildDSTree()
            
            issat = 0
            if t != False:
              issat = 1
              
            geneSpeciesMapping = getGeneSpeciesMapping(genes, speciesTree, ":")
              
            graph = ConstraintGraph(genes, orthologs, paralogs, speciesTree, geneSpeciesMapping)
            t = graph.buildDSTree()
            
            iscons = 0
            if t != False:
              iscons = 1
              
            print str(prob_orth) + "," + str(prob_para) + "," + str(nbGenes) + "," + str(nbSpecies) + "," + str(issat) + "," + str(iscons)
'''

    
