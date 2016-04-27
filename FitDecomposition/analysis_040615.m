%% analyze Nichola's updated data
%original version: 04/06/2015
%updated: 11/19/2015
%this script is similar to analysis_031715
%used slightly different data (from Nick)
%function: data processing; find complete subgraphs (findinter,findinter_index)

clear;
clc;

folder='./analysis_040615/';
load mutfitness

%% data processing
%load original data
% load([folder 'Mutfit_updated.mat']);
% %filter
% %1)filter low input (NaN in I20fit, Nicholas already filtered the data for input <10)
% index_input=find(~isnan(Mutfit.I20fit)); 
% %2)filter stop codon
% index_nostop=find(cellfun(@isempty, strfind(Mutfit.mut,'_')));
% %entries left
% index=intersect(index_input,index_nostop);
%filtered data
% mutfitness=Mutfit(index,:);

% %set lethal mutations to 10^-4 (this was only used in Fourier analysis, log(fitness))
% index_lethal=find(mutfitness.I20fit==0);
% mutfitness.I20fit(index_lethal)=10^-4;

% save
% save('mutfitness.mat','mutfitness');

%% find WT index
WT_seq='VDGV';
index_WT =find(ismember(mutfitness.mut,WT_seq));

%% find mutants with Hamming Distance=4, and their intermediates
index_HD{4}=find(mutfitness.HD==4);
%group into 1,2,3 mutants
index_HD{3}=find(mutfitness.HD==3);
index_HD{2}=find(mutfitness.HD==2);
index_HD{1}=find(mutfitness.HD==1);

%% find genotypes of intermediates
total=length(index_HD{4});
tic;
index_start=1;
index_end=total;

graph_genotype=cell(total,1);
for i=index_start:index_end
    graph_genotype{i}=findinter(WT_seq, char(mutfitness.mut(index_HD{4}(i))));
end
toc;

%% Hamming distance for all genotypes in the subgraph (L=4, A=2)
L=4;
graph_total=2^L;
genotype_str=dec2base(0:graph_total-1,2);
genotype_num=str2num_sequence(genotype_str);
Hdistance=pdist2(genotype_num(1,:),genotype_num,'Hamming')*L;

%% find index/fitness of intermediate genotypes
tic;
[graph_index, complete]=findinter_index(mutfitness.mut, graph_genotype, index_HD, index_WT,Hdistance);
toc;
%Elapsed time is 1392.061772 seconds.

%% save
%information of complete subgraphs
% save('graph_genotype.mat','graph_genotype');
% save('graph_index.mat','graph_index');
% save('graph_other.mat','complete','Hdistance','index_WT','index_HD','index_lethal');

