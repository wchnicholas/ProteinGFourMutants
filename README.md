This README describes the scrips used for the analyses in: 
[Adaptation in protein fitness landscapes is facilitated by indirect paths](https://elifesciences.org/content/5/e16965/article-metrics)

### FILES
* doc/Epistasis: Epistasis value from [Olson et al. 2014](http://www.cell.com/current-biology/abstract/S0960-9822(14)01268-8)
* doc/SMutList:  Single mutant read count from [Olson et al. 2014](http://www.cell.com/current-biology/abstract/S0960-9822(14)01268-8)
* doc/DMutList:  Double mutant read count from [Olson et al. 2014](http://www.cell.com/current-biology/abstract/S0960-9822(14)01268-8)
* doc/AAmass: Mass and volume for each amino acid
* Fasta/Barcode.fa: Barcodes used for sequencing.
  * Input: cDNA input library
  * IGG10: low concentration of IgG is used as target (unpublished results)
  * IGG20: medium concentration of IgG is used as target
  * IGG90: high concentration of IgG is used as target (unpublished results)
* Fasta/SeqInfo.fa: reference sequence information
* result/regression\_missing: inferred fitness for missing variants in log scale
* count/: This folder contains the number of occurrence for each variant in the sequencing data
  * Input.count: cDNA input library
  * IGG10.count: low concentration of IgG is used as target (unpublished results)
  * IGG20.count: medium concentration of IgG is used as target
  * IGG90.count: high concentration of IgG is used as target (unpublished results)

### MAPPING PIPELINE
* Mapper1\_NW.py: Extract nucleotide information from read and demultiplexing. Files for forward reads must contain '\_R1\_'. Files for reverse reads must contain '\_R2\_'.
  * Input file: 
    * fastq/\*.fastq
    * Fasta/SeqInfo.fa
    * Fasta/Barcode.fa
  * Output file:
    * paired/Input 
    * paired/IGG10 
    * paired/IGG20 
    * paired/IGG90 
* Mapper2\_NW.py: Translating nucleotide sequence to protein sequence
  * Input file: 
    * paired/Input 
    * paired/IGG10 
    * paired/IGG20 
    * paired/IGG90 
  * Output file: 
    * pep/Input.pep
    * pep/IGG10.pep
    * pep/IGG20.pep
    * pep/IGG90.pep
* Mapper3\_NW.sh: Counting the number for each unique protein sequence
  * Input file: 
    * pep/Input.pep
    * pep/IGG10.pep
    * pep/IGG20.pep
    * pep/IGG90.pep
  * Output file: 
    * count/Input.count
    * count/IGG10.count
    * count/IGG20.count
    * count/IGG90.count
* Mapper4\_NW.py: Combine count information from different samples into one file, also compute hamming distance
  * Input file: 
    * count/Input.count
    * count/IGG10.count
    * count/IGG20.count
    * count/IGG90.count
  * Output file: 
    * result/WTcount
    * result/Mutcount
* Mapper5\_NW.py: Compute fitness information, and also add the column for Anders single/double mutant fitness data, and Kd calculation
  * Input file: 
    * result/WTcount
    * result/Mutcount
  * Output file: 
    * result/Mutfit

### ANALYSIS
* Analysis1a\_NW.py: Compute episasis for all partition combinations for HD = 4 (Unused)
  * Input file:  result/Mutfit
  * Output file: result/HD4EpiIGG\*

* Analysis1b\_NW.py: Compute episasis for using S1\*S2\*S3\*S4 for all mutations. (Unused)
  * Input file:  result/Mutfit
  * Output file: result/AllEpi

* Analysis2\_NW.py: Framework for graph analysis and finding local maximum, finding accessible paths, etc.

* Analysis3a\_NW.py: Pathway analysis, e.g. count number of shortest pathways from WT to beneficial mutations
  * Input file:  result/Mutfit
  * Output file: analysis/ShortPaths4muts\*

* Analysis3b\_NW.py: Pathway analysis, e.g. count number of shortest pathways from deleterious variants to WT. Restricted to subgraph with one-localmax (WT being the sole localmax). First step has to be > 0.01
  * Input file:  result/Mutfit
  * Output file: analysis/ShortMonoPaths4ToWT\*

* Analysis4\_NW.py: Find fitness effect of a particular mutant in different backgrounds
  * Input file:  result/Mutfit
  * Output file:
    * analysis/MutDiffBG\*
    * analysis/FitDiffBG 

* Analysis5\_NW.py: Adapted from Analysis2.py. Print pathway parameters for each variant

* Analysis6\_NW.py: Extract shortest pathways from data, converting to xdot/graphviz, can tune the start node and end node

* Analysis7\_NW.py: Search for epistasis under different backgrounds. Search for positive in one and negative in other.
  * Input file:  result/Mutfit
  * Output file: analysis/EpiDiffBGI\*

* Analysis8\_TMP.py: Calculate epistasis under different models, including epistatic accumulation model

* Analysis9\_NW.py: Search for three way interactions (Unused)
  * Input file:  result/Mutfit
  * Output file: analysis/ThreeWayEpi

* Analysis10a\_NW.py: Calculate Gini Index, Bias of Pathway (From WT to beneficial mutations)
  * Input file:  result/Mutfit
  * Output file: analysis/PathwayParamResult

* Analysis10b\_NW.py: Calculate Gini Index, Bias of Pathway (From fit < 1 mutants to WT)
  * Input file:  result/Mutfit
  * Output file: analysis/PathwayParamResultToWT

* Analysis11\_NW.py: Global graph analysis, global max search, pathway to max search
  * Input file:  result/Mutfit
  * Output file:
    * analysis/LocalMaxMuts\*
    * analysis/LocalMaxClimb\*
    * analysis/LocalMaxCompile\*

* Analysis11sim\_NW.py: A high-throughput computing version of Analysis11\_NW.py

* Analysis12\_NW.py: Plot the distance relationship between local max, fitness cutoff > 1, edge = HD 2
  * Input file:  analysis/LocalMaxCompile
  * Output file:
    * xdot/LocalMax.png
    * xdot/LocalScale.png

* Analysis12b\_NW.py: Plot the distance relationship between local max from both inferred 2nd order landscape and original landscape
  * Input file:
    * analysis/LocalMaxCompile\_greedy
    * analysis/LocalMaxCompile\_greedy\_pair
  * Output file:
    * xdot/LocalMax\_combine.png
    * xdot/LocalScale\_combine.png

* Analysis13\_NW.py: Find minimum pathlength (monotonic) to those local max with a fitness of > 1
  * Input file:
    * result/Mutfit
    * result/regression\_missing
    * analysis/LocalMaxCompile\_greedy
  * Outputfile: analysis/LocalMaxPathLen

* Analysis14\_NW.py: Compute the path length information (i.e. avg, min, max) from the simulation output
  * Input file:
    * simulations/\*/LocalMaxClimb\*
    * analysis/LocalMaxDes\_weight\_\*
  * Output file:  analysis/LocalMaxDist\_weight\_\*

* Analysis15\*\_NW.py: Calculate the evolution potential of each variant (how many variants it can reach by upward path)
  * Input file:
    * result/Mutfit
    * result/regression_missing
  * Output file:
    * analysis/LocalMaxEvolvePotSampled (Analysis15\_NW.py)
    * analysis/LocalMaxEvolvePotR\*      (Analysis15R\_NW.py)
    * analysis/LocalMaxEvolvePotWT      (Analysis15WT\_NW.py)

* Analysis15WTnuc\_NW.py: Calculate the evolution potential of WT with the contraints from standard genetic code
  * Input file: 
    * result/Mutfit
    * result/regression_missing
  * Output file: 
    * analysis/AAtransitionmatrix
    * analysis/LocalMaxEvolvePotWTnuc

* Analysis16\_NW.py: Compute the neighbor correlation of the entire landscape (Unused)

* Analysis17.py: Search for rerouting paths

* Analysis18\_NW.py: Compile Mass and Fitness data
  * Input file: 
    * result/Mutfit
    * doc/AAmass
  * Output file: analysis/BulkvsFit

* Analysis19\_NW.py: Mutational Step Classifier
  * Input file: analysis/LocalMaxClimb\_greedy or simulations/\*/LocalMaxClimb\_\*
  * Output file: It produces a standard output, which is copied into analysis/StepsInfo

* Analysis20\_NW.py: Bootstrapping the average hamming distances among randomly sampled variants

* Analysis21\_NW.py: Compute the fitness difference between neighboring variants (Unused)
  * Input file: 
    * result/Mutfit
    * result/regression\_missing

* Analysis22.py: Compile the endpoint reproducibility for each variant (Unused)
  * Input file:
    * simulations/\*/LocalMaxClimb\_\*
    * analysis/LocalMaxDes\_\*
  * Output file: analysis/RepTraj\_\*

* Analysis23.py: Count accessible peaks for a given node and count how many nodes can reach a given peak
  * Input file:  analysis/LocalMaxPathLen

* Analysis24.py: Customize search for epistasis pair of interest
  * Input file: result/Mutfit 
  * Output file: analysis/AllPairwiseEpi

* EvolvePotFromWT.sh: Compute evolution potential and accessbility to beneficial mutations from WT

### FITNESS DECOMPOSITION
scripts in FitDecomposition/ are written by Lei Dai. They are for fitness decomposition analysis using fourier transform. 

* FitDecomp1.py: Compute fitness information for each subgraphs
  * Input file:  analysis/FitnessDecompose
  * Output file: analysis/FitnessDecomposeFit

### PARSING ADAPTIVE PATHWAY SIMULATIONS
* CompileProbDest.py: Compile simulation results in to a file, count the evolution endpoint, pathway reproducibility and diversity
  * Input file:
    * result/Mutfit
    * result/regression\_missing
    * simulations/\*/LocalMaxClimb\_\*
  * Output file: analysis/LocalMaxDes\_\*

### PLOTTING SCRIPTS
* Plot1\_NW.R: Plot correlation between conditions and also against Anders single/double mutant fitness data
  * Input file: result/Mutfit
  * Output file: graph/G4Cor\_\*.png

* Plot2a\_NW.R: Plot DFE as boxplot and hist for different IGG conditions and different HD groups 
  * Input file: result/Mutfit
  * Output file:
    * graph/DFEboxIGG\*.png
    * graph/DFEhistIGG\*.png

* Plot2b\_NW.R: Plot Distribution of Epistasis, fitness and expected fitness as histogram for different IGG conditions and different HD groups 
  * Input file:  result/AllEpi
  * Output file:
    * graph/DEEIGG\*.png
    * graph/DEEhistIGG\*.png
    * graph/ExpfitboxIGG\*.png 

* Plot3\_NW.R: Plot the maximum predicted fitness across predictions from different types of partition vs actual fitness for HD = 4.
  * Input file: result/HD4EpiIGG20

* Plot4\_NW.R: Plot the maximum and minimum fitness across different genetic backgrounds for individual single substituition
  * Input file:
    * analysis/MutDiffBGI20fit
    * analysis/FitDiffBG
  * Output file:
    * graph/DiffBGIGG20.png
    * graph/DFEDiffBG.png

* Plot5\_NW.R: Plot the maximum and minimum epistasis, and also give some stats about the frequency of such phenomenon
  * Input file:  analysis/EpiDiffBGI\*
  * Output file: graph/ContextEpi\*.png

* Plot6\_NW.R: Plot the correlation result from fitness decomposition
  * Input file:  analysis/FitnessDecompose
  * Output file:
    * graph/FitnessDecompHist.png
    * graph/FitnessDecomp.png

* Plot7\_NW.R: Plot the results from pathway analysis, e.g. number of stucking pathways, number of permissive pathways
  * Input file:
    * analysis/ShortPaths4mutsCutoff10fold
    * analysis/ShortPaths4mutsCutoff1fold
  * Output file:
    * graph/PathwayDistribute.png
    * graph/PathwayDistribute\_stuck.png
    * graph/PathwayDistribute\_mono.png

* Plot8.R: Unfinished/Unused

* Plot9.R: Unfinished/Unused

* Plot10.R: Plot the relationship between basin size and fitness of the local maximums and also plot the number of step to local maximums from each variant
  * Input file: 
    * analysis/LocalMaxCompile\_\*
    * analysis/LocalMaxClimb\_\*
  * Output file:
    * graph/LocalMaxBasinvsFit\_\*.png
    * graph/LocalMaxPathtoMax\_greedy.png

* Plot11.R: Plot reproducibility, entropy, evolutionary potential
  * Input file:
    * analysis/LocalMaxPathLen
    * analysis/LocalMaxMuts
    * analysis/LocalMaxDes\_random
  * Output file:
    * graph/LocalMax15accessbox.png
    * graph/LocalMax15accesshist.png
    * graph/LocalMaxDesReprodAll\_\*.png
    * graph/LocalMaxPathReprodAll\_\*.png
    * graph/LocalMaxDesReprodBenMut\_\*.png
    * graph/LocalMaxPathReprodBenMut\_\*.png
    * graph/LocalMaxDes\_lowpeak.png

* Plot12.R: Plot the fraction of reachable beneficial variants. (Unused)
  * Input file:
    * analysis/LocalMaxEvolveAll
    * analysis/LocalMaxEvolveAll\_pair
  * Output file: graph/LocalMaxEvolvePot

* Plot13.R: Plot entropy along mutational trajectories
  * Input file: analysis/RepTraj\_weight
  * Output file: graph/TrackingEn\_weight.png

* Plot14.R: Plot the heatmap for epistasis of a given substitution pair under different genetic backgrounds

* ManFiguring.R: Plotting figures for manuscript

* Heatmapping1.py: Format the data into a heatmap format, for EpiRange and EpiSD
  * Input file:  analysis/EpiDiffBGI20fit
  * Output file:
    * analysis/HeatMapEpiRange
    * analysis/HeatMapEpiSD

* Heatmapping2.py: Format the data into a heatmap format, for a given mutation pair under all possible genetic backgrounds
  * Input file:  result/Mutfit
  * Output file: analysis/Heatmap\_\*

* Heatmapping3.py: Clustering mutant (Unused)
  * Input file: analysis/LocalMaxDes\_\*
  * Outputfile: analysis/LocalMaxDesCluster\_\*

### MISC SCRIPTING
* BasicInfo1.R: For extraction of basic statistics from the data (e.g. Total coverage, Maximum fitness, etc.)
  * Input file:  result/Mutfit
  * Output file: BasicInfo.ods and FitnessDiffHD.ods
* BasicInfo2.R: Computing the correlation between different predictors
  * Input file:  result/HD4EpiIGG\*
  * Output file  NONE
* OPTM.R: For optimization of Bmax
* Grepping all path lengths into a tmp file
  * ``awk {'print $1, $2'} simulations/LocalMaxClimb_weight* | grep -v 'steps' > tmp/SimPaths``
* Printing beneficial mutants according to hamming distance
  * awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 1) print $1}' result/Mutfit > BenMutHD1
  * awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 2) print $1}' result/Mutfit > BenMutHD2
  * awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 3) print $1}' result/Mutfit > BenMutHD3
  * awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 4) print $1}' result/Mutfit > BenMutHD4

