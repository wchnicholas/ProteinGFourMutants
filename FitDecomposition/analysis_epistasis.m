%last edit: 11/16/2015
clear;
clc;
%load original data
%fitness
load analysis_040615/mutfitness
%genotype
%do not use genotype_int from regression! for regression, I changed the mapping such that WT 'VDGV'=[20 20 20 20]
load bypass/genotype_int

%%
% for i=1:height(mutfitness)
% genotype_int(i,:)=aa2int(mutfitness.mut{i});
% end
% save('./bypass/genotype_int','genotype_int');

%% randomly sample a genotype, two single mutations and the double mutation
N=size(genotype_int,1); %number of observations
L=4; %genome length
sample=1e5;
% sample=100;
count=0;

tic;
while count<sample
    %sample genotype
    index_sample=randsample(N,1);
    genotype_sample_int= genotype_int(index_sample,:);
    genotype_sample=int2aa(genotype_sample_int);
    %sample mutations
    mutation_pos = datasample(1:L,2,'Replace',false); %site # of two mutations, random draw without replacement
    %residue 1
    residue_1=genotype_int(index_sample,mutation_pos(1));
    temp=1:1:20;
    mutationspec_1=temp(temp~=residue_1);
    mutation_1=randsample(mutationspec_1,1);
    %residue 2
    residue_2=genotype_int(index_sample,mutation_pos(2));
    temp=1:1:20;
    mutationspec_2=temp(temp~=residue_2);
    mutation_2=randsample(mutationspec_2,1);
    
    %single mutants
    genotype_single1_int=genotype_sample_int;
    genotype_single1_int(mutation_pos(1))=mutation_1;
    genotype_single1=int2aa(genotype_single1_int);
    genotype_single2_int=genotype_sample_int;
    genotype_single2_int(mutation_pos(2))=mutation_2;
    genotype_single2=int2aa(genotype_single2_int);
    %double mutant
    genotype_double_int=genotype_sample_int;
    genotype_double_int(mutation_pos(1))=mutation_1;
    genotype_double_int(mutation_pos(2))=mutation_2;
    genotype_double=int2aa(genotype_double_int);
    %find genotypes in the list
    index_single1 =find(ismember(mutfitness.mut,genotype_single1));
    index_single2 =find(ismember(mutfitness.mut,genotype_single2));
    index_double =find(ismember(mutfitness.mut,genotype_double));
    if ~isempty(index_single1) && ~isempty(index_single2) && ~isempty(index_double)
        index_temp(1)=index_sample;
        index_temp(2)=index_single1;
        index_temp(3)=index_single2;
        index_temp(4)=index_double;
        fitness_temp=mutfitness.I20fit(index_temp);
        genotype_temp=double(genotype_int(index_temp,:));
        
        %exclude samples with >=2 lethal mutants: they are set to equal values, this causes problem in epistasis classification
        exclude=0;
        %two lethal mutants, but not adajcent
        lethal=find(fitness_temp==1e-4);
        if length(lethal)==2 && pdist2(genotype_temp(lethal(1),:),genotype_temp(lethal(2),:),'hamming')==1/L
            exclude=1;
        else if length(lethal)>2
                exclude=1;
            end
        end
        if ~exclude
            count=count+1;
            index_all{count}=index_temp;
        end
    end
end
toc;
%~30s for 1000 samples

%% classify the type of epistasis
epistasis_count=zeros(3,1);
count_reciprocal=0;
tic;
for i=1:sample
    genotype_quad=double(genotype_int(index_all{i},:));
    fitness_quad=mutfitness.I20fit(index_all{i});
    %classify pairwise epistasis: epistasis_classify.m
    epistasis_type=epistasis_classify(fitness_quad,genotype_quad);
    epistasis_count(epistasis_type)= epistasis_count(epistasis_type)+1;
    
    %find reciprocal epistasis
    if epistasis_type==3
        count_reciprocal=count_reciprocal+1;
        index_reciprocal{count_reciprocal}=index_all{i};
    end

end
toc;

%% save
% save('./bypass/index_reciprocal_1e5','index_all','index_reciprocal');

%% calculate the fraction of each type
epistasis_freq=epistasis_count/sample;
%plot
bar(epistasis_freq);
ylabel('Fraction of pairwise epistasis');
set(gca,'XTickLabel',{'Magnitude', 'Sign', 'Reciprocal sign'})

