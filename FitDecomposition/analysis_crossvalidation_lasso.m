%% 10-fold cross validation for a range of lambda: determine the optimal penalty
%http://www.mathworks.com/help/stats/lasso.html
clear;
clc;

%inputs
%predictor_all_order3: pairwise interactions + 39-41-54 three way interaction
%total:9102 parameters
%fitness_all: fitness values of all genotypes, input<10 are not included
%total:119884 observations

%%  
%load original data
%fitness
load analysis_040615/mutfitness
%genotype_int: only used for regression, shifted VDGV to 20 20 20 20
load ./regression/genotype_int_shifted

%predictor matrix (generate each time, because it's too large to store)
%20-25 GB memory needed during the matrix generation (could be optimized by combining predictor_all with 3rd order terms)
%parameters
categ=[1 2 3 4];
catlevels=[20 20 20 20];
%genotype
genotype_all=genotype_int;
%fitness
fitness_all=log(mutfitness.I20fit);

%exclude lethal mutation in fitting
index_nonlethal=find(mutfitness.I20fit>10^-4);
genotype_all=genotype_int(index_nonlethal,:);
fitness_all=fitness_all(index_nonlethal);

predictor_all=x2fx(genotype_all,'interaction',categ, catlevels);  %designer matrix to predictor matrix

%include third-order interaction terms in predictor matrix
%limit to position 39,41,54 (position 1,3,4)
genotype_reduced=double(genotype_all(:,[1;3;4]));
index_order3=genotype_reduced(:,3)+(genotype_reduced(:,2)-1)*19+(genotype_reduced(:,1)-1)*19^2;
predictor_order3=zeros(length(index_order3),19^3);
for i=1:length(index_order3)
    if max(genotype_reduced(i,:))<20 %mutations at three residues
         predictor_order3(i,index_order3(i))=1;
    end
end

%combine with predictor matrix with 1st and 2nd order terms (x2fx)
num_sample=size(predictor_all,1);
num_predictor=size(predictor_all,2)+size(predictor_order3,2);
predictor_all_order3=zeros(num_sample,num_predictor);
predictor_all_order3=[predictor_all, predictor_order3];

%clear memory
predictor_all=[];
predictor_order3=[];

%% lasso, penalty=1e-4,3e-3,1e-3,3e-2,1e-2
penalty_CV=[1e-4 3e-3 1e-3 3e-2 1e-2]; %I meant 3e-4...

%the CV option in lasso funtion turns out to be infeasible because memory crashed.
%tic;
%[beta_nonlethal_lasso_cv, FitInfo_nonlethal_lasso_cv]=lasso(predictor_all_order3,fitness_all,'lambda',penalty_cv,'CV',10);
%toc;

%% generate training and test set 
%cross validation, 10-fold
K=10;
% N=height(mutfitness);
N=length(fitness_all);
indices = crossvalind('Kfold', N, K);
%
save('./regression/lasso_nonlethal_CVpartition','indices');

%%
for j=1:length(penalty_CV)
    tic;
    for i = 1:K
        test = (indices == i);
        train = ~test;          
        %lasso
        [beta_nonlethal_lasso_CV{j,i}, FitInfo_nonlethal_lasso_CV{j,i}]=lasso(predictor_all_order3(train,:),fitness_all(train),'lambda',penalty_CV(j));
    end
    toc;
    %save after finishing each penalty parameter
    save('./regression/lasso_nonlethal_CV','beta_nonlethal_lasso_CV','FitInfo_nonlethal_lasso_CV');
end

%% this run was finished over Thanksgiving
% Elapsed time is 185738.386457 seconds.
% Elapsed time is 33614.068721 seconds.
% Elapsed time is 87199.535870 seconds.
% Elapsed time is 18012.678566 seconds.
% Elapsed time is 20111.153074 seconds.

