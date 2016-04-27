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
%20-25 GB memory needed during the matrix generation
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

%%
%combine with predictor matrix with 1st and 2nd order terms (x2fx)
predictor_order3=[predictor_all, predictor_order3];

%clear memory
predictor_all=[];

%%
% load ./regression/beta_nonlethal_order3_CV
load ./regression/crossvalidation_nopenalty %11/15/2015 CV run (need to run again)
load ./regression/lasso_nonlethal_CV

%% MSE: regression, 10-fold cross validation (analysis_MSE.m)
K=10;
for i=1:K
    MSE_CV(i)=mean((fitness_test_regression{i}-fitness_test{i}).^2);
end
mean_MSE=mean(MSE_CV);
std_MSE=std(MSE_CV);

%% plot MSE vs penalty: lasso
figure(1);
num_penalty=size(FitInfo_nonlethal_lasso_CV,1);
num_fold=size(FitInfo_nonlethal_lasso_CV,2);
for i=1:num_penalty
    lambda(i)=FitInfo_nonlethal_lasso_CV{i,1}.Lambda;
    MSE_temp=[];
    for j=1:num_fold
        MSE_temp=[MSE_temp FitInfo_nonlethal_lasso_CV{i,j}.MSE];
    end
    mean_MSE_lasso(i)=mean(MSE_temp);
    std_MSE_lasso(i)=std(MSE_temp);
end
errorbar(lambda,mean_MSE_lasso,std_MSE_lasso,'o');

line([10^-4 10^-1],[mean_MSE mean_MSE],'color','k');
line([10^-4 10^-1],[mean_MSE+std_MSE mean_MSE+std_MSE],'color','r');
line([10^-4 10^-1],[mean_MSE-std_MSE mean_MSE-std_MSE],'color','r');

set(gca,'xscale','log','xlim',[10^-4 10^-1]);
xlabel('Penalty parameter (\lambda)');
ylabel('10-fold cross-validation Mean Squared Error');
% set(gca,'fontsize',15);

%% plot nonzero elements
figure(2);
nnz_temp=[];
for i=1:K
    nnz_temp=[nnz_temp nnz(beta_nonlethal_order3_CV{i})];
end
mean_nnz=mean(nnz_temp);
std_nnz=std(nnz_temp);

for i=1:num_penalty
    nnz_temp=[];
     for j=1:num_fold
         nnz_temp=[nnz_temp nnz(beta_nonlethal_lasso_CV{i,j})];
     end
     mean_nnz_lasso(i)=mean(nnz_temp);
     std_nnz_lasso(i)=std(nnz_temp);
end
errorbar(lambda,mean_nnz_lasso,std_nnz_lasso,'o');
    
    
line([10^-4 10^-1],[mean_nnz mean_nnz],'color','k');
line([10^-4 10^-1],[mean_nnz+std_nnz mean_nnz+std_nnz],'color','r');
line([10^-4 10^-1],[mean_nnz-std_nnz mean_nnz-std_nnz],'color','r');
set(gca,'xscale','log','xlim',[10^-4 10^-1]);
xlabel('Penalty parameter (\lambda)');
ylabel('Number of nonzero parameters');
% set(gca,'fontsize',15);
%% choose penalty=1e-4
