%% fit  linear regression model to fitness of measured variants
%goal: impute fitness of missing variants
%1) include single site and pairwise interaction terms
%2) include third-order interaction terms: a) residue 39,41,54; b) all
%last update: 11/14/2015

clear;
%load original data
%fitness
load analysis_040615/mutfitness
%genotype_int: only used for regression, shifted VDGV to 20 20 20 20
load ./regression/genotype_int_shifted

%load previous fit
% load analysis_040615/regression_tenfold

%% transform sequence to integers
% tic;
% for i=1:height(mutfitness)
% genotype_int(i,:)=aa2int(mutfitness.mut{i}); %in regression 0000 would be VVVV 
% 
% %quick fix: 0000 would be the real WT, VDGV, modified to expand around VDGV
% if mutfitness.mut{i}(2)=='D' %swap D<->V
%     genotype_int(i,2)=20;
% else if mutfitness.mut{i}(2)=='V'
%         genotype_int(i,2)=4;
%     end
% end
% 
% if mutfitness.mut{i}(3)=='G' %swap D<->V
%     genotype_int(i,3)=20;
% else if mutfitness.mut{i}(3)=='V'
%         genotype_int(i,3)=8;
%     end
% end
% 
% end
% toc;
%Elapsed time is 61.753581 seconds.
% save('./regression/genotype_int_shifted','genotype_int');

%% fit with all data
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

%% predictor matrix
predictor_all=x2fx(genotype_all,'interaction',categ, catlevels);  %designer matrix to predictor matrix

%% check predictor matrix: find nonzero terms
% %check is successful
% %WT: VDGV
% %single: ADGV; double: AAGV; triple; AAAV; quadraple: AAAA
% test_seq='VDGV';
% % test_seq='AAAA';
% index_test =find(ismember(mutfitness.mut,test_seq));
% genotype_int(index_test,:)
% find(predictor_all(index_test,:))

%% include third-order interaction terms in predictor matrix
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

%%
% %% regression
% tic;
% % beta_all=regress(fitness_all,predictor_all);
% beta_all_order3=regress(fitness_all,predictor_all_order3);
% toc;
% % Warning: X is rank deficient to within machine precision. 
% % > In regress at 84 
% % Elapsed time is 3790.038200 seconds.
% 
% %% lasso
% % penalty=logspace(-6,-4,3);
% penalty=10^-2;
% for i=1:length(penalty)
%     tic;
%     [beta_all_lasso2{i}, FitInfo2{i}]=lasso(predictor_all,fitness_all,'lambda',penalty(i));
%     toc;
% end

%% plot fit
%all data
% fitness_all_regression=(predictor_all*beta_all);
% fitness_all_regression=(predictor_all_order3*beta_all_order3);
fitness_all_regression=(predictor_all_order3*beta_nonlethal_order3);
correlation_all=corr(fitness_all,fitness_all_regression,'type','Pearson');

figure(11);
set(gcf, 'Visible', 'on');
plot(fitness_all, fitness_all_regression,'+');
legend(strcat('Pearson=',num2str(correlation_all)));
title('All data');
xlabel('True fitness');
ylabel('Predicted fitness');

%% MSE of linear regression
% fitness_all_regression=(predictor_all*beta_all);
MSE=mean((fitness_all_regression-fitness_all).^2);

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
                
%% predict fitness of all/missing genotypes
predictor_comb=x2fx(genotype_comb_int,'interaction',categ, catlevels);  
fitness_comb=(predictor_comb*beta_all);

%% find missing genotypes
[genotype_missing_int, index_missing]=setdiff(genotype_comb_int,genotype_all, 'rows');
genotype_missing=cell(length(index_missing),1);
for i=1:length(index_missing)
  genotype_missing{i}=int2aa(genotype_missing_int(i,:));
end

%%
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


%% combine fitness data: measured value of available gentoypes + predicted value of missing genotypes 
%04/20/15
%sort index of available genotypes
index_available=zeros(size(genotype_int,1),1);
for i=1:size(genotype_int,1)
    temp=0;
    for j=1:4
        temp=temp+(20^(4-j))*double(genotype_int(i,j)-1); %in this case, genotype_int is from 1 to 20
    end
    index_available(i)=temp+1;
end

index_missing=zeros(size(genotype_missing_int,1),1);

%sort index of missing genotypes
for i=1:size(genotype_missing_int,1)
    temp=0;
    for j=1:4
        temp=temp+(20^(4-j))*double(genotype_missing_int(i,j)-1); %in this case, genotype_int is from 1 to 20
    end
    index_missing(i)=temp+1;
end

fitness_combined=zeros(20^4,1);
fitness_combined(index_available)=fitness_all; %data (log transformed)
fitness_combined(index_missing)=fitness_comb(index_missing); %predicted values

%%
save('fitness_combined_VVVV','fitness_combined');

%% output: send to Nick, 04/06/15
T_regression_all=table(genotype_comb,  fitness_comb);
T_regression_missing=table(genotype_missing, fitness_comb(index_missing));

% writetable(T_regression_all,'regression_all.dat');
% writetable(T_regression_missing,'regression_missing.dat');
 
%% save
% save('regression_tenfold.mat','genotype_int','indices','beta','correlation_train','correlation_test','beta_all');

%% explore: rank vs fitness

fitness_comb_sorted=sort(fitness_comb,'descend');
fitness_union=[log(mutfitness.I20fit); fitness_comb(index_missing)];
fitness_sorted=sort(fitness_union,'descend');

%%
figure(100);
plot(fitness_sorted/log(10),'b');
hold on;
plot(fitness_comb_sorted/log(10),'r');
set(gca,'xscale','log');

%% analyze coefficients
hist(beta_all);
xlabel('coefficients');

