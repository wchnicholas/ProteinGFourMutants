%% plot: bar plot for 3 types of epistasis (Figure 1C)
%11/24/2015
clear;
clc;

%% around WT
load ./bypass/index_reciprocal_WT
% calculate the fraction of each type
epistasis_count_WT=zeros(3,1);
for i=1:length(epistasis_type_WT)
 epistasis_count_WT(epistasis_type_WT(i))= epistasis_count_WT(epistasis_type_WT(i))+1;
end
epistasis_freq_WT=epistasis_count_WT/length(epistasis_type_WT);

%% randomly sampled pairs in the sequence space
load analysis_040615/mutfitness
load bypass/genotype_int
load ./bypass/index_reciprocal_1e5
%classify the type of epistasis
epistasis_count=zeros(3,1);
tic;
for i=1:sample
    genotype_quad=double(genotype_int(index_all{i},:));
    fitness_quad=mutfitness.I20fit(index_all{i});
    %classify pairwise epistasis: epistasis_classify.m
    epistasis_type=epistasis_classify(fitness_quad,genotype_quad);
    epistasis_count(epistasis_type)= epistasis_count(epistasis_type)+1;
end
toc;
%calculate the fraction of each type
sample=1e5;
epistasis_freq=epistasis_count/sample;

%% save
% save('./bypass/epistasis_freq','epistasis_freq','epistasis_freq_WT');

%% plot
epistasis_plot=[epistasis_freq_WT';epistasis_freq']';
b=bar(epistasis_plot);
b(1).LineWidth = 2;
b(1).EdgeColor = 'b';
b(1).FaceColor = 'b';
b(2).LineWidth = 2;
b(2).EdgeColor = [1 0.5 0.5];
b(2).FaceColor = [1 0.5 0.5];
ylabel('Fraction');
legend('Around WT','Entire sequence space');
set(gca,'ylim',[0 1]);
box off;
set(gca,'XTickLabel',{'Magnitude', 'Sign', 'Reciprocal sign'})
set(gca,'fontsize',15');
% legend('Magnitude epistasis','Sign epistasis','Reciprocal sign epistasis');