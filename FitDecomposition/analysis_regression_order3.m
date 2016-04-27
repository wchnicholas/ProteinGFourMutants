%% set up predictor matrix: analysis_regression.m

%% fit filtered variants 
%regression
tic;
beta_nonlethal_order3=regress(fitness_all,predictor_all_order3);
toc;
save('./regression/beta_nonlethal_order3','beta_nonlethal_order3');

%lasso  
penalty=logspace(-4,-1,4);
for i=1:length(penalty)
    tic;
    [beta_nonlethal_lasso{i}, FitInfo_nonlethal{i}]=lasso(predictor_all_order3,fitness_all,'lambda',penalty(i));
    toc;
end
save('./regression/lasso_nonlethal_order3','beta_nonlethal_lasso','FitInfo_nonlethal');

%% 
%regression
% Warning: X is rank deficient to within machine precision. 
% > In regress at 84 
% Elapsed time is 2977.734989 seconds.

%lasso, penalty=10^-4,-3,-2,-1
% Elapsed time is 21083.236002 seconds. 
% Elapsed time is 9050.756250 seconds.
% Elapsed time is 1939.734821 seconds.
% Elapsed time is 1401.651119 seconds.


