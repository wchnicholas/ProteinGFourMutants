%
clear;
%load original data
%fitness
load analysis_040615/mutfitness
%genotype_int: only used for regression, shifted VDGV to 20 20 20 20
load ./regression/genotype_int_shifted

%% load previous fit
% load ./regression/beta_nonlethal_order3

load ./regression/lasso_nonlethal_order3 %from analysis_regression_order3.m
%correlation coefficients: penalty=1e-4
beta_lasso=beta_nonlethal_lasso{1};
%add intercept: corrected on 02/11/2016
intercept=FitInfo_nonlethal{1}.Intercept;

%% predictor matrix (generate each time, because it's too large to store)
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
fitness_all_regression=(predictor_all_order3*beta_nonlethal_order3);
%lasso
fitness_all_lasso=(predictor_all_order3*beta_lasso)+intercept;
correlation_all=corr(fitness_all,fitness_all_regression,'type','Pearson');

%%
figure(1);
subplot(1,2,1);
set(gcf, 'Visible', 'on');
plot(fitness_all, fitness_all_regression,'+');
legend(strcat('Pearson=',num2str(correlation_all)));
set(gca,'fontsize',15);
title('All data');
xlabel('Measured fitness');
ylabel('Predicted fitness');
%draw x=y line
hold on;
temp_plot=-10:0.01:5;
plot(temp_plot,temp_plot,'-r');

subplot(1,2,2);
%plot histogram of residuals
hist(fitness_all_regression-fitness_all);

%% lasso prediction
figure(2);
% subplot(1,2,1);
set(gcf, 'Visible', 'on');
plot(fitness_all, fitness_all_lasso,'+');
legend(strcat('Pearson=',num2str(correlation_all)));
set(gca,'fontsize',15);
% title('All data');
xlabel('Measured fitness');
ylabel('Predicted fitness');
%draw x=y line
hold on;
temp_plot=-10:0.01:4;
plot(temp_plot,temp_plot,'-r');
set(gca,'xlim',[-10 4],'ylim',[-10 4]);

%%
% subplot(1,2,2);
%plot histogram of residuals
hist(fitness_all_lasso-fitness_all_regression);
set(gca,'fontsize',15);
xlabel('Residuals');

%% load previous fit
load ./regression/lasso_nonlethal_order3
load ./regression/lasso_nonlethal_order3_run2

%% lasso: plot correlation
for i=1:4
fitness_all_lasso=(predictor_all_order3*beta_nonlethal_lasso{i});
correlation_lasso=corr(fitness_all,fitness_all_lasso,'type','Pearson');

subplot(2,2,i);
set(gcf, 'Visible', 'on');
plot(fitness_all, fitness_all_lasso,'+');
legend(strcat('Pearson=',num2str(correlation_lasso)));
set(gca,'fontsize',15);
title(strcat('Penalty=',num2str(FitInfo_nonlethal{i}.Lambda)));
xlabel('Measured fitness');
ylabel('Predicted fitness');
end

%% MSE of linear regression
fitness_all_regression=(predictor_all_order3*beta_nonlethal_order3);
MSE=mean((fitness_all_regression-fitness_all).^2);

%% plot MSE vs penalty
for i=1:4
    plot(FitInfo_nonlethal{i}.Lambda,FitInfo_nonlethal{i}.MSE,'o');
    hold on;
end
line([10^-4 10^-1],[MSE MSE],'color','k');
set(gca,'xscale','log','xlim',[10^-4 10^-1]);
set(gca,'fontsize',15);
xlabel('Penalty');
ylabel('MSE');