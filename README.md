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