%% analyze pairwise epistasis around WT
%updated: 11/24/2015
clear;
clc;
%fitness
load analysis_040615/mutfitness
%genotype_int
load bypass/genotype_int

%% find reciprocal sign epistasis, WT as starting point
WT_seq='VDGV';
WT_seq_int=aa2int(WT_seq);
index_sample=find(ismember(mutfitness.mut,WT_seq));
genotype_sample_int=WT_seq_int;

count_pair=0;
count=0;
epistasis_type_WT=[];
%enumerate pairs of mutations
mutation_pos = combnk(1:4,2); %site # of two mutations

tic;
for i=1:size(mutation_pos,1);
    %residue 1
    residue_1=genotype_int(index_sample,mutation_pos(i,1));
    temp=1:1:20;
    mutationspec_1=temp(temp~=residue_1);
    %residue 2
    residue_2=genotype_int(index_sample,mutation_pos(i,2));
    temp=1:1:20;
    mutationspec_2=temp(temp~=residue_2);
    
    for j=1:length(mutationspec_1) %enumerate all amino acids
        for k=1:length(mutationspec_2)
            mutation_1=mutationspec_1(j);
            mutation_2=mutationspec_2(k);
            %single mutants
            genotype_single1_int=genotype_sample_int;
            genotype_single1_int(mutation_pos(i,1))=mutation_1; %bug fixed. 11/18/2015
            genotype_single1=int2aa(genotype_single1_int);
            genotype_single2_int=genotype_sample_int;
            genotype_single2_int(mutation_pos(i,2))=mutation_2;
            genotype_single2=int2aa(genotype_single2_int);
            %double mutant
            genotype_double_int=genotype_sample_int;
            genotype_double_int(mutation_pos(i,1))=mutation_1;
            genotype_double_int(mutation_pos(i,2))=mutation_2;
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
                
                count_pair=count_pair+1;
                %classify pairwise epistasis: epistasis_classify.m
                epistasis_type_WT(count_pair)=epistasis_classify(fitness_temp,genotype_temp);
  
                %identify reciprocal epistasis (WT as the starting point)
                if fitness_temp(1)>fitness_temp(2) && fitness_temp(1)>fitness_temp(3) && fitness_temp(4)>fitness_temp(2) && fitness_temp(4)>fitness_temp(3)
                    count=count+1;
                    index_reciprocal_WT{count}=index_temp;                        
                end                
            end
            
        end
    end
end
toc;

%% save
% save('./bypass/index_reciprocal_WT','epistasis_type_WT','index_reciprocal_WT');

%% plot: fraction of pairwise epistasis types around WT
% calculate the fraction of each type
epistasis_count_WT=zeros(3,1);
for i=1:length(epistasis_type_WT)
 epistasis_count_WT(epistasis_type_WT(i))= epistasis_count_WT(epistasis_type_WT(i))+1;
end
epistasis_freq_WT=epistasis_count_WT/length(epistasis_type_WT);
bar(epistasis_freq_WT);
ylabel('Fraction');
set(gca,'XTickLabel',{'magnitude epistasis', 'sign epistasis', 'reciprocal sign epistasis'})


%% check
% index_single=find(mutfitness.HD==1);
% subplot(1,2,1);
% hist(mutfitness.I20fit(index_single));
% index_double=find(mutfitness.HD==2);
% subplot(1,2,2);
% hist(mutfitness.I20fit(index_double));

%% check
% test_seq='VDFV';
% test_seq='WDGV';
% test_seq='WDFV';
% index_test=find(ismember(mutfitness.mut,test_seq));
% test_fitness=mutfitness.I20fit(index_test)