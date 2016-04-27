%12/5/2015
%export data for Nick's simulations
%reference: analysis_regression_compare (expand around VDGV, used shifted int2aa mapping)

clear
clc;
load ./regression/lasso_nonlethal_order3 %from analysis_regression_order3.m

%correlation coefficients: penalty=1e-4
beta_lasso=beta_nonlethal_lasso{1};
%add intercept: corrected on 02/11/2016
intercept=FitInfo_nonlethal{1}.Intercept;

%% enumerate all genotypes
genotype_comb_int=zeros(20^4,4);
genotype_comb=cell(20^4,1);
count=1;
for i=1:20
    for j=1:20
        for k=1:20
            for l=1:20
                genotype_comb_int(count,:)=[i j k l];
                genotype_comb{count}=int2aa([i j k l]);              
                count=count+1;
            end
        end
    end
end

%% see analysis_regression: genotype_int_shifted (this section is not necessary)
%position2: swap D<->V; position 3: swap G<->V
for i=1:20^4
    if genotype_comb_int(i,2)== 20 %swap D<->V
        genotype_comb{i}(2)='D';
    else if genotype_comb_int(i,2)== 4
            genotype_comb{i}(2)='V';
        end
    end
    
    if genotype_comb_int(i,3)== 20 %swap G<->V
        genotype_comb{i}(3)='G';
    else if genotype_comb_int(i,3)==8
            genotype_comb{i}(3)='V';
        end
    end
  
end
%% predictor matrix
%parameters
categ=[1 2 3 4];
catlevels=[20 20 20 20];
predictor_comb=x2fx(genotype_comb_int,'interaction',categ, catlevels);

%% include third-order interaction terms in predictor matrix
%limit to position 39,41,54 (position 1,3,4)
genotype_reduced=double(genotype_comb_int(:,[1;3;4]));
index_order3=genotype_reduced(:,3)+(genotype_reduced(:,2)-1)*19+(genotype_reduced(:,1)-1)*19^2;
predictor_order3=zeros(length(index_order3),19^3);
for i=1:length(index_order3)
    if max(genotype_reduced(i,:))<20 %mutations at three residues
         predictor_order3(i,index_order3(i))=1;
    end
end

%combine predictor matrix 
predictor_order3=[predictor_comb, predictor_order3];
%clear memory
predictor_comb=[];

%% predict fitness of all/missing genotypes
fitness_comb=(predictor_order3*beta_lasso)+intercept; %corrected, 02/11/2016
%save
% save('./regression/fitness_comb_lasso1e-4_corrected','fitness_comb');

%% find missing genotypes
%genotype observed (not missing)
load ./regression/genotype_int_shifted
genotype_all=genotype_int;
[genotype_missing_int, index_missing]=setdiff(genotype_comb_int,genotype_all, 'rows');
genotype_missing=cell(length(index_missing),1);
for i=1:length(index_missing)
  genotype_missing{i}=int2aa(genotype_missing_int(i,:));
end

%% transform missing genotypes from int to aa alphabet: this is important ('VDGV'=[20 20 20 20])
for i=1:length(index_missing)
  genotype_missing{i}=int2aa(genotype_missing_int(i,:));
  
if genotype_missing_int(i,2)== 20 %swap D<->V
    genotype_missing{i}(2)='D';
else if genotype_missing_int(i,2)== 4
        genotype_missing{i}(2)='V';
    end
end

if genotype_missing_int(i,3)== 20 %swap G<->V
    genotype_missing{i}(3)='G';
else if genotype_missing_int(i,3)==8
        genotype_missing{i}(3)='V';
    end
end
  
end

%% output: corrected aa2int mapping
% T_regression_all=table(genotype_comb,  fitness_comb);
% writetable(T_regression_all,'regression_all_WT.dat');

T_regression_missing=table(genotype_missing, fitness_comb(index_missing));
writetable(T_regression_missing,'./regression/lasso1e-4_missing_corrected.dat');

%% check: histogram of fitness values of missing genotypes 
histogram(fitness_comb(index_missing));
