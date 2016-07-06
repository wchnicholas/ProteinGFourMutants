#MAPPING PIPELINE
* Mapper1_NW.py: Extract nucleotide information from read and demultiplexing
* Mapper2_NW.py: Translating nucleotide sequence to protein sequence
* Mapper3_NW.sh: Counting the number for each unique protein sequence
* Mapper4_NW.py: Combine count information from different samples into one file, also compute hamming distance
* Mapper5_NW.py: Compute fitness information, and also add the column for Anders single/double mutant fitness data, and Kd calculation

###########################################PLOTTING SCRIPTS##################################################
Plot1_NW.R: Plot correlation between conditions and also against Anders single/double mutant fitness data
  #Input file:  result/Mutfit
  #Output file: graph/G4Cor_*.png

Plot2a_NW.R: Plot DFE as boxplot and hist for different IGG conditions and different HD groups 
  #Input file:  result/Mutfit
  #Output file: graph/DFEboxIGG*.png
  #             graph/DFEhistIGG*.png

Plot2b_NW.R: Plot Distribution of Epistasis, fitness and expected fitness as histogram for different IGG conditions and different HD groups 
  #Input file:  result/AllEpi
  #Output file: graph/DEEIGG*.png
  #             graph/DEEhistIGG*.png
  #             graph/ExpfitboxIGG*.png 

Plot3_NW.R: Plot the maximum predicted fitness across predictions from different types of partition vs actual fitness for HD = 4.
  #Input file:  result/HD4EpiIGG20

Plot4_NW.R: Plot the maximum and minimum fitness across different genetic backgrounds for individual single substituition
  #Input file:  analysis/MutDiffBGI20fit
  #             analysis/FitDiffBG
  #Output file: graph/DiffBGIGG20.png
  #             graph/DFEDiffBG.png

Plot5_NW.R: Plot the maximum and minimum epistasis, and also give some stats about the frequency of such phenomenon
  #Input file:  analysis/EpiDiffBGI*
  #Output file: graph/ContextEpi*.png

Plot6_NW.R: Plot the correlation result from fitness decomposition (Lei Dai result)
  #Input file:  analysis/FitnessDecompose
  #Output file: graph/FitnessDecompHist.png
  #             graph/FitnessDecomp.png

Plot7_NW.R: Plot the results from pathway analysis, e.g. number of stucking pathways, number of permissive pathways
  #Input file:  analysis/ShortPaths4mutsCutoff10fold
  #             analysis/ShortPaths4mutsCutoff1fold
  #Output file: graph/PathwayDistribute.png
  #             graph/PathwayDistribute_stuck.png
  #             graph/PathwayDistribute_mono.png

Plot8.R: Unfinished

Plot9.R: Unfinished

Plot10.R: 1. Plot the relationship between basin size and fitness of the local maximums.
          2. Plot the number of step to local maximums from each variant
  #Input file: analysis/LocalMaxCompile_*
  #            analysis/LocalMaxClimb_*
  #Output file: graph/LocalMaxBasinvsFit_*.png
  #             graph/LocalMaxPathtoMax_greedy.png

Plot11.R: Plot reproducibility, entropy, evolutionary potential
  #Input file: analysis/LocalMaxPathLen
  #            analysis/LocalMaxMuts
  #             analysis/LocalMaxDes_random
  #Output file: graph/LocalMax15accessbox.png
  #             graph/LocalMax15accesshist.png
  #             graph/LocalMaxDesReprodAll_*.png
  #             graph/LocalMaxPathReprodAll_*.png
  #             graph/LocalMaxDesReprodBenMut_*.png
  #             graph/LocalMaxPathReprodBenMut_*.png
  #             graph/LocalMaxDes_lowpeak.png

Plot12.R: Plot the fraction of reachable beneficial variants. 
  #Input file: analysis/LocalMaxEvolveAll
  #            analysis/LocalMaxEvolveAll_pair
  #Output file: graph/LocalMaxEvolvePot

Plot13.R: Plot entropy along mutational trajectories
  #Input file: analysis/RepTraj_weight
  #Output file: graph/TrackingEn_weight.png

Plot14.R: Plot the heatmap for epistasis of a given substitution pair under different genetic backgrounds

ManFiguring.R: Plotting figures for manuscript

Heatmapping1.py: Format the data into a heatmap format, for EpiRange and EpiSD
  #Inputfile:  analysis/EpiDiffBGI20fit
  #Outputfile: analysis/HeatMapEpiRange
  #            analysis/HeatMapEpiSD

Heatmapping2.py: Format the data into a heatmap format, for a given mutation pair under all possible genetic backgrounds
  #Inputfile:  result/Mutfit
  #Outputfile: analysis/Heatmap_*

Heatmapping3.py: Clustering mutant (Unused)
  #Input file: analysis/LocalMaxDes_*
  #Outputfile: analysis/LocalMaxDesCluster_*

############################################ANALYSIS SCRIPTS####################################################
Analysis1a_NW.py: Compute episasis for all partition combinations for HD = 4
  #Input file:  result/Mutfit
  #Output file: result/HD4EpiIGG*

Analysis1b_NW.py: Compute episasis for using S1*S2*S3*S4 for all mutations, THIS SCRIPT HAS A GOOD EPISTASIS CALCULATION
  #Input file:  result/Mutfit
  #Output file: result/AllEpi

Analysis2_NW.py: Framework for graph analysis and finding local maximum, finding accessible paths, etc.

Analysis3a_NW.py: Pathway analysis, e.g. count number of shortest pathways from WT to beneficial mutations
  #Input file:  result/Mutfit
  #Output file: analysis/ShortPaths4muts*

Analysis3b_NW.py: Pathway analysis, e.g. count number of shortest pathways from deleterious variants to WT
                  Restricted to subgraph with one-localmax (WT being the sole localmax)
                  First step has to be > 0.01
  #Input file:  result/Mutfit
  #Output file: analysis/ShortMonoPaths4ToWT*

Analysis4_NW.py: Find fitness effect of a particular mutant in different backgrounds
  #Input file:  result/Mutfit
  #Output file: analysis/MutDiffBG*
  #             analysis/FitDiffBG 

Analysis5_NW.py: Adapted from Analysis2.py. Print pathway parameters for each variant

Analysis6_NW.py: Extract shortest pathways from data, converting to xdot/graphviz, can tune the start node and end node

Analysis7_NW.py: Search for epistasis under different backgrounds. Search for positive in one and negative in other.
  #Input file:  result/Mutfit
  #Output file: analysis/EpiDiffBGI*

Analysis8_TMP.py: Calculate epistasis under different models, including epistatic accumulation model

Analysis9_NW.py: Search for three way interactions (Unused)
  #Input file:  result/Mutfit
  #Output file: analysis/ThreeWayEpi

Analysis10a_NW.py: Calculate Gini Index, Bias of Pathway (From WT to beneficial mutations)
  #Input file:  result/Mutfit
  #Output file: analysis/PathwayParamResult

Analysis10b_NW.py: Calculate Gini Index, Bias of Pathway (From fit < 1 mutants to WT)
  #Input file:  result/Mutfit
  #Output file: analysis/PathwayParamResultToWT

Analysis11_NW.py: Global graph analysis, global max search, pathway to max search
  #Input file:  result/Mutfit
  #Output file: analysis/LocalMaxMuts
  #             analysis/LocalMaxClimb
  #             analysis/LocalMaxCompile

Analysis11sim_NW.py: A high-throughput computing version of Analysis11_NW.py

Analysis12_NW.py: Plot the distance relationship between local max, fitness cutoff > 1, edge = HD 2
  #Input file:  analysis/LocalMaxCompile
  #Output file: xdot/LocalMax.png
  #             xdot/LocalScale.png

Analysis12b_NW.py: Plot the distance relationship between local max from both inferred 2nd order landscape and original landscape
  #Input file:  analysis/LocalMaxCompile_greedy
  #             analysis/LocalMaxCompile_greedy_pair
  #Output file: xdot/LocalMax_combine.png
  #             xdot/LocalScale_combine.png

Analysis13_NW.py: Find minimum pathlength (monotonic) to those local max with a fitness of > 1
  #Input file:  result/Mutfit
  #             result/regression_missing
  #             analysis/LocalMaxCompile_greedy
  #Outputfile:  analysis/LocalMaxPathLen

Analysis14_NW.py: Compute the path length information (i.e. avg, min, max) from the simulation output
  #Input file:   simulations/*/LocalMaxClimb*
  #              analysis/LocalMaxDes_weight_pair
  #Output file:  analysis/LocalMaxDist_weight_pair

Analysis15*_NW.py: Caculate the evolution potential of each variant (how many variants it can reach by upward path)
  #Input file:   result/Mutfit
  #              result/regression_missing
  #Output file:  analysis/LocalMaxEvolvePotSampled (Analysis15_NW.py)
  #              analysis/LocalMaxEvolvePotR*      (Analysis15R_NW.py)
  #              analysis/LocalMaxEvolvePotWT      (Analysis15WT_NW.py)

Analysis16_NW.py: Compute the neighbor correlation of the entire landscape

Analysis17.py: Search for rerouting paths

Analysis18_NW.py: Compile Mass and Fitness data
  #Input file:   result/Mutfit'
  #              doc/AAmass
  #Output file:  analysis/BulkvsFit'

Analysis19_NW.py: Mutational Step Classifier
  #Input file: analysis/LocalMaxClimb_greedy or simulations/*/LocalMaxClimb_*
  #Output file: It produces a standard output, which is copied into analysis/StepsInfo

Analysis20_NW.py: Bootstrapping the average hamming distances among randomly sampled variants

Analysis21_NW.py: Compute the fitness difference between neighboring variants
  #Input file:  result/Mutfit
  #             result/regression_missing
  #Output file: analysis/FitStepDist

Analysis22.py: Compile the endpoint reproducibility for each variant
  #Input file: simulations/*/LocalMaxClimb_*
  #            analysis/LocalMaxDes_*
  #Output file: analysis/RepTraj_*

Analysis23.py: Count accessible peaks for a given node and count how many nodes can reach a given peak
  #Input file:  analysis/LocalMaxPathLen
  #Output file: Doc/NodesAccessPeaks.ods

Analysis24.py: Customize search for epistasis pair of interest
  #Input file: result/Mutfit 
  #Output file: analysis/AllPairwiseEpi

EvolvePotFromWT.sh: Compute evolution potential and accessbility to beneficial mutations from WT
  #Output file: ../Doc/EvolvePotFromWT.ods

#################FITNESS DECOMPOSITION###################
scripts in FitDecomposition are written by Lei Dai. They are for fitness decomposition analysis using fourier transform. 

FitDecomp1.py: Compute fitness information for each subgraphs
  #Input file:  analysis/FitnessDecompose
  #Output file: analysis/FitnessDecomposeFit

#################PATHWAY SIMULATIONS#####################

CompileProbDest.py: Compile simulation results in to a file, count the evolution endpoint, pathway reproducibility and diversity
  #Input file:  result/Mutfit
  #             result/regression_missing
  #             simulations/*/LocalMaxClimb_*
  #Output file: analysis/LocalMaxDes_weight

#################To grep all pathway lengths into a tmp file##############
awk {'print $1, $2'} simulations/LocalMaxClimb_weight* | grep -v 'steps' > tmp/SimPaths

for m in `awk {'if ($4 == -1) print $1'} analysis/LocalMaxEvolvePotWT_pair`; do grep ^$m analysis/LocalMaxEvolvePotWT; done > tmp/pathofinterest


awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 1) print $1}' result/Mutfit > BenMutHD1
awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 2) print $1}' result/Mutfit > BenMutHD2
awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 3) print $1}' result/Mutfit > BenMutHD3
awk '{if ($11 > 1 && $11 != "NA" && $1 != "mut" && $2 <= 4) print $1}' result/Mutfit > BenMutHD4

#MISC SCRIPTING
* BasicInfo1.R: For extraction of basic statistics from the data (e.g. Total coverage, Maximum fitness, etc.)
  * Input file:  result/Mutfit
  * Output file: BasicInfo.ods and FitnessDiffHD.ods
* BasicInfo2.R: Computing the correlation between different predictors
  * Input file:  result/HD4EpiIGG\*
  * Output file  NONE
* OPTM.R: For optimization of Bmax
* Grepping all path lengths into a tmp file

	awk {'print $1, $2'} simulations/LocalMaxClimb_weight* | grep -v 'steps' > tmp/SimPaths

#FILES
result/Epistasis: Epistasis value from Anders et al. 2014
