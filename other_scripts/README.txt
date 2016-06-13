The scripts in this directory were used to obtain the results discussed in the paper.
They are shamelessly undocumented and worse, they contain hard-coded paths and parameters.
The author blames time advancing too quickly for the unforgivable sins.


If you are interested in the parameters given to proetinortho, they are hardcoded in jsonExtractor.py


Nonetheless, anyone interested in reproducing/extending the results can modify and use them
(though it would be best to contact lafonman@iro.umontreal.ca - I am willing to offer my collaboration
to navigate this pipeline).


Here is a brief description of what the scripts accomplish : 

--- jsonExtractor.py ---

Fetches a list of Ensembl gene tree IDs from a hard-coded gene tree file and, 
for each of them, obtains the nucleic sequences (using Ensembl json rest api) of each gene.
These sequences are sent to proteinortho using 5 hard-coded parameter sets from loose to strict,
generating 5 graph files in the current directory.




--- results_maker.py ---
Fetches a list of Ensembl gene tree IDs from a hard-coded gene tree file and,
for each of them, checks whether the five proteinortho graph files are present.
If so, runs constraint_checker.py on each of them individually, and then on every pair
of graph files to check for satisfiability and consistency with a hard-coded species tree.
The results are sent to a .satcons file, one per tree.



--- ensembl_rels_fetcher.py ---
For every .satcons file that's found in a hard-coded directory, manage to fetch the gene tree
used to construct it and extract all orthology/paralogy relationships depicted by this tree
(according to LCA reconciliation).  The relationships go to a .ensrels file.


--- results_parser.py ---
Reads every .satcons file and their corresponding .ensrels file, and extracts some numbers out of it.
Results are divided in 2 categories : single graph by type, double graph by type-type.
For single graphs, results look like 
--------ultrastrict-----------
NBISSAT=55
NBISCONS=18

In this example, the number of satisfiable, then consistent constraint graphs under ultrastrict orthology inference.

Double graphs go like
--------loose:strict-----------
AVG_PCT_BADORTH_CONS=0.151962486212
AVG_PCT_BADORTH_DUB_CONS=0.0509476376333
AVG_PCT_BADORTH_DUB_SAT=0.052774927318
AVG_PCT_BADORTH_SAT=0.179473345164
AVG_PCT_BADPARA_CONS=0.171776698103
AVG_PCT_BADPARA_SAT=0.209718153862
AVG_PCT_REL_CONS=0.506454685349
AVG_PCT_REL_SAT=0.527168727474
NBISCONS=149
NBISSAT=254
NB_CONS_WHEN_NO_SINGLE_CONS=109
NB_SAT_WHEN_NO_SINGLE_SAT=142
SUM_PCT_BADORTH_CONS=22.6424104456
SUM_PCT_BADORTH_DUB_CONS=7.59119800736
SUM_PCT_BADORTH_DUB_SAT=13.4048315388
SUM_PCT_BADORTH_SAT=45.5862296716
SUM_PCT_BADPARA_CONS=25.5947280174
SUM_PCT_BADPARA_SAT=53.268411081
SUM_PCT_REL_CONS=75.4617481169
SUM_PCT_REL_SAT=133.900856778

Let's explain those numbers : 
AVG_PCT_BADORTH_CONS : Of all consistent trees, the average of (orthologies in both graphs and in Ensembl)/(orthologies in both graphs)
AVG_PCT_BADORTH_DUB_CONS= Of all consistent trees, the average of (orthologies in both graphs and in Ensembl + not dubious in Ensembl)/(orthologies in both graphs)0.0509476376333
AVG_PCT_BADORTH_DUB_SAT=Of all satisfiable trees, the average of (orthologies in both graphs and in Ensembl + not dubious in Ensembl)/(orthologies in both graphs)0.0509476376333
AVG_PCT_BADORTH_SAT=Of all satisfiable trees, the average of (orthologies in both graphs and in Ensembl)/(orthologies in both graphs)0.0509476376333
AVG_PCT_BADPARA_CONS=Same as AVG_PCT_BADORTH_CONS, but for paralogies
AVG_PCT_BADPARA_SAT=Same as AVG_PCT_BADORTH_SAT, but for paralogies
AVG_PCT_REL_CONS=The average number of DECIDED (ie not undecided) relationships when consistent
AVG_PCT_REL_SAT=The average number of DECIDED (ie not undecided) relationships when satisfiable
NBISCONS=Number of consistent graphs in this mode
NBISSAT=Number of satisfiable graphs in this mode
NB_CONS_WHEN_NO_SINGLE_CONS=Number of consistent graphs in this mode, when no single graph was consistent (ie we saved the situation)
NB_SAT_WHEN_NO_SINGLE_SAT=Number of satisfiable graphs in this mode, when no single graph was satisfiable (ie we saved the situation)

And the sums are just temporary variable needed to compute the averages.