%% cross validation: linear regression without penalty
%updated: 11/16/2015

%% cross validation, 10-fold
K=10;
% N=height(mutfitness);
N=length(fitness_all);
indices = crossvalind('Kfold', N, K);

%% parameters (used for generating predictor matrix, pairwise interaction)
% categ=[1 2 3 4];
% catlevels=[20 20 20 20];

%% fit
tic;
for i = 1:K
test = (indices == i);
train = ~test;

%genotype
genotype_train=genotype_int(train,:);
%fitness
fitness_train=fitness_all(train);
%predictor matrix
% predictor_train=x2fx(genotype_train,'interaction',categ, catlevels);  %designer matrix to predictor matrix
predictor_train=predictor_all_order3(train,:);

%regression
beta_nonlethal_order3_CV{i}=regress(fitness_train,predictor_train);
end
toc;
save('./regression/beta_nonlethal_order3_CV','beta_nonlethal_order3_CV');


%% archive: previous run with pairwise interactions
% genotype_test=genotype_int(test,:);
% predictor_test=x2fx(genotype_test,'interaction',categ, catlevels);  
%10 regression: Elapsed time is 1684.564294 seconds. 


%% calculate correlation
beta=beta_nonlethal_order3_CV;
tic;
for i=1:K
    %true fitness of test data set (validation)
    test = (indices == i);
    fitness_test{i}=fitness_all(test);
    predictor_test=predictor_all_order3(test,:);
    
    % evaluate: correlation between true and predicted fitness
    fitness_test_regression{i}=(predictor_test*beta{i});
    correlation_test(i)=corr(fitness_test{i},fitness_test_regression{i},'type','Pearson');
end
toc;

%% plot 
for i=1:K
    subplot(2,5,i);
    set(gcf, 'Visible', 'on');
%     figure(i);
%     subplot(1,2,1);
%     plot(fitness_train, fitness_train_regression,'+');
%     legend(strcat('Pearson=',num2str(correlation_train(i))));
%     title('Train data');
%     xlabel('True fitness');
%     ylabel('Predicted fitness');
%     subplot(1,2,2);
    plot(fitness_test{i}, fitness_test_regression{i},'+');
    legend(strcat('Pearson=',num2str(correlation_test(i))));
    title('Test data');
    xlabel('True fitness');
    ylabel('Predicted fitness');
    
%     filename=strcat(pwd,'\figure\regression\correlation_set',num2str(i));
%     print('-djpeg ','-cmyk',filename);

end

%% save analysis, 11/15/2015 run
% save('./regression/crossvalidation_nopenalty',...
%     'fitness_all','genotype_int','indices','beta_nonlethal_order3_CV','fitness_test','fitness_test_regression','correlation_test');
