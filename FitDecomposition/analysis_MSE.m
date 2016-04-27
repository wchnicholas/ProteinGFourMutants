clear;
clc;
load crossvalidation_nopenalty
load beta_nonlethal_order3
load beta_nonlethal_order3_CV
load lasso_nonlethal_order3
load lasso_nonlethal_order3_run2

%% MSE: regression, 10-fold cross validation
K=10;
for i=1:K
    MSE_CV(i)=mean((fitness_test_regression{i}-fitness_test{i}).^2);
end
mean_MSE=mean(MSE_CV);
std_MSE=std(MSE_CV);

%% plot MSE vs penalty: lasso
figure(1);
num_penalty=size(FitInfo_nonlethal,2);
for i=1:num_penalty
    plot(FitInfo_nonlethal{i}.Lambda,FitInfo_nonlethal{i}.MSE,'o');
    hold on;
end
num_penalty=size(FitInfo_nonlethal_run2,2);
for i=1:num_penalty
    plot(FitInfo_nonlethal_run2{i}.Lambda,FitInfo_nonlethal_run2{i}.MSE,'o');
    hold on;
end

line([10^-4 10^-1],[mean_MSE mean_MSE],'color','k');
line([10^-4 10^-1],[mean_MSE+std_MSE mean_MSE+std_MSE],'color','r');

set(gca,'xscale','log','xlim',[10^-4 10^-1]);
xlabel('Penalty');
ylabel('MSE');

%% plot nonzero elements
figure(2);
num_penalty=size(FitInfo_nonlethal,2);
for i=1:num_penalty
    plot(FitInfo_nonlethal{i}.Lambda,nnz(beta_nonlethal_lasso{i}),'o');
    hold on;
end
num_penalty=size(FitInfo_nonlethal_run2,2);
for i=1:num_penalty
    plot(FitInfo_nonlethal_run2{i}.Lambda,nnz(beta_nonlethal_lasso_run2{i}),'o');
    hold on;
end

line([10^-4 10^-1],[nnz(beta_nonlethal_order3) nnz(beta_nonlethal_order3)],'color','k');
line([10^-4 10^-1],[length(beta_nonlethal_order3) length(beta_nonlethal_order3)],'color','r');
set(gca,'xscale','log','xlim',[10^-4 10^-1]);
xlabel('Penalty');
ylabel('# of nonzero parameters');