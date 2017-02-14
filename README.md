## BeWith v0.1

**BeWith** is a general clustering framework for identifying modules with different combinations of mutation and interaction patterns. On a high level, given a set of genes and two types of edge scoring functions (within and between scores), **BeWith** aims to find clusters of genes so that genes within a cluster maximize the “within” scores while gene pairs spanning two different clusters maximize the “between” scores.
There are **three** different settings of the BeWith framework: 
+ **bemewithfun**: ensuring mutual exclusivity of mutations between different modules and functional similarity of genes within modules.
+ **bemewithco**: ensuring mutual exclusivity between modules and co-occurrence of mutations in genes within modules.
+ **becowithmefun**: ensuring co-occurrence between modules while enforcing mutual exclusivity and functional interactions within modules. 


### Requirements

+ Linux/Mac OS/Windows
+ CPLEX version 12+
+ Apache Maven version 3+ (optional for recompilation of the code)

**Please make sure CPLEX (cplex or cplex.exe) is callable from the folder that you intend to run BeWith**

### Compilation

Prepackaged jar file is available in the target folder. If you decide to recompile, please make sure you have Apache Maven installed. 

Download the source code and change into package folder:
```
git clone https://github.com/phuongdao1/bewith.git
cd bewith
```
Produce the package jar file using Maven:
```
mvn package
```
If packaging process is successful, the jar file *BeWithFramework.jar* should be produced inside folder *target*. 

### Usage

Usage: java -jar BeWithFramework.jar -m METHOD -k NUMBER_OF_CLUSTERS -o OUTPUT_PREFIX [-NETWORK_TYPE1 INPUT_FILE1] [-NETWORK_TYPE2 INPUT_FILE2]...

* METHOD can be one of the following:    
  + *bemewithfun*, to ensure mutual exclusivity of mutations between different modules and functional similarity of genes within modules
  + *bemewithco*, to ensure mutual exclusivity between modules and co-occurrence of mutations in genes within modules
  + *becowithmefun*, to ensure co-occurrence between modules while enforcing mutual exclusivity and functional interactions within modules

* NUMBER OF CLUSTERS should be at least 2 for bemewithfun and bemewithco, currently always set at 2 for *becowithfun*

* OUTPUT_PREFIX is the prefix of the result and intermediate files

* NETWORK_TYPE can be of the following:        
  + *fun*, for a functional or protein interaction network. The input file for the functional/protein network should consist of three space/tab-seperated columns. The gene names are in the first two columns and the score/weights (between 0 and 1) of the interactions are in the last columns. Note that **BeWith** assumes the interactions are undirectional. Here is an example:
  
  ```
  M6PR    PLIN3   0.903
  ESRRA   SIRT1   0.946
  ESRRA   NRF1    0.956
  ESRRA   NRBP1   0.911
  ESRRA   NR1D1   0.909
  ...
  ```
  + *me*, for a network built from mutually exclusive mutationns. The input file for this network also consists of three space/tab-seperated columns. The gene names are in the first two columns and the p-values indicating the significance of the relations are in the last columns. Followed is an example:
  ```
  LRBA  PIK3CA  0.0455
  UBR5  PIK3CA  0.0434
  MST1P9  PIK3CA  0.037012
  NBPF1 PIK3CA  0.00821
  NCOR2 MUC16 0.0105
  ...
  ```
  + *co*, for a network built from co-occurring mutations. The input file for the mutational co-occurrence network also consists of three space/tab-seperated columns. The gene names are in the first two columns and the p-values indicating the significance of the relations are in the last columns. Here is an example:
  ```
  COL4A6  COL6A3  0.007878034268846874
  OR2T34  PLEC  0.0055867103553033015
  WWP2  THSD7A  0.004376064336983943
  ERBB3 COL1A2  0.0043760643369844246
  ANKRD12 MYLK  0.008597917892609233
  ...
  ```

### A Walkthrough Example

Here we will go through an example of running BeWith framework. The provided dataset is TCGA BRCA cancer data and is processed as described in the [Bewith manuscript](#citation). Pairs of genes with significant p-values for mutually exclusive mutations computed from [WeSME](#other-references) are stored *data/BRCA_me_WESME.txt*. Similarly, pairs of genes with significant p-values for co-occurring mutations by hypergeometric tests are stored in *data/BRCA_co.txt*. Gene pairs with significant p-values for co-occurring mutations computed from [WeSCO](#other-references) are stored in *data/BRCA_co_WESCO.txt. The functional protein STRING network file *HumanStringNet.txt* is also included in folder *data*.

We utilize **BeWith** to find co-occurrence modules which are mutually exclusive to each other. In other words, genes within each modules have co-occurring mutation and pairs of genes belong to different modules likely have mutually exclusive mutations. In what follows, we utilize **BeWith** to find five such modules using *bemewithco* option: 

```
java -jar BeWithFramework.jar -m bemewithco -k 5 -o BRCA_bemewithco -me data/BRCA_me_WESME.txt -co data/BRCA_co.txt
```

If there are no errors, the summary of the results should be in the file *BRCA_bemewithco_5_modules* and looks like this:

```
Score = 37.650293769959475

Module 1: TP53 DCC
Module 2: MUC16 TTN GON4L
Module 3: PCDH19 GATA3
Module 4: PIK3CA MAP2K4
Module 5: NCOA3 NCOR2 MEF2A

Within-module edges:
PIK3CA,MAP2K4,co,1.0590367
NCOA3,MEF2A,co,3.0
NCOR2,NCOA3,co,3.0
NCOR2,MEF2A,co,2.8316467
MUC16,TTN,co,1.4000814
MUC16,GON4L,co,1.1947811
TTN,GON4L,co,1.7089661
TP53,DCC,co,1.3005306
GATA3,PCDH19,co,1.0006765

Between-module edges:
PIK3CA,TP53,me,3.0
PIK3CA,PCDH19,me,0.68125516
PIK3CA,DCC,me,0.87527853
PIK3CA,MEF2A,me,0.68125516
PIK3CA,GON4L,me,0.7306503
PIK3CA,GATA3,me,2.2101083
NCOR2,MUC16,me,0.98940533
MUC16,TP53,me,0.7592787
MUC16,GATA3,me,1.1302139
TTN,NCOA3,me,1.0228788
TP53,NCOA3,me,1.1858056
TP53,MAP2K4,me,2.0378604
TP53,GATA3,me,3.0
DCC,TTN,me,0.76362175
GATA3,TTN,me,2.0869627
```

The summary report starts with the objective function value of the solution of the ILP model. Then comes the list of modules together with genes within the modules. We list pairs of genes with significant mutally exclusive/co-occurring mutations or with functional interactions (applicable in other two settings *bemewithfun* and *becowithmefun*) in two sections: pairs of genes in the same module *Within-module edges* and pairs of genes in different modules *Between-module edges*. *me*/*co* denotes mutually exclusive/co-occurring mutations. The scores of the gene pairs are computed as minus of log to the base of 100 of the p-values. Score of gene pairs with p-values less than 1/10^6 are set to 3.0. *pp* denotes a functional interaction or a protein interaction gene pair.     

### Citation

If you use BeWith in your work, please cite:

Phuong Dao\*, Yoo-Ah Kim\*, Damian Wojtowicz, Sanna Madan, Roded Sharan, and Teresa M. Przytycka. BeWith: A Between-Within Method to Discover Relationships between  Cancer Modules via  Integrated Analysis of Mutual Exclusivity, Co-occurrence and Functional Interactions. To appear.

\* equal contribution


### Other References

Yoo-Ah Kim  Sanna Madan  Teresa M. Przytycka. WeSME: Uncovering Mutual Exclusivity of Cancer Drivers and Beyond. Bioinformatics (2016) btw242. DOI: https://doi.org/10.1093/bioinformatics/btw242.
