%% analyze Fourier spectrum, Pearson correlation
%last update: 11/25/2015
clear;
clc;

folder='./analysis_040615/';
%complete graphs
load([folder 'mutfitness']);
load([folder 'graph_index']);
load([folder 'graph_genotype']);
load([folder 'graph_other']);
load([folder 'graph_fitness']);

%saved analysis
load([folder 'graph_pearson']);
load([folder 'graph_spectrum']);
% load graph_secondorder

%find complete subgraphs
index_complete=find(complete==1);

%% spectrum CDF
% total_plot=length(index_complete);
L=4;
total_plot=100;
figure(100);
for i=1:total_plot
n=index_complete(i);
plot(graph_spectrum{n},'--','color',[0.5 0.5 0.5]);
hold on;
end

set(gca,'xlim',[1 L]);
xlabel('Order of Fourier decomposition');
ylabel('Spectrum CDF');

%% second order expansion: compare specturm and Pearson correlation
spectrum_plot=[];
pearson_plot=[];
for i=1:length(index_complete)
    n=index_complete(i);
   spectrum_plot(n)= graph_spectrum{n}(2);
    pearson_plot(n)=graph_pearson{n}(2);
end

%% histogram/CDF of spectrum and Pearson correlation
figure(100);
subplot(2,2,1);
hist(spectrum_plot(index_complete));
set(gca,'xlim',[0 1]);
title('second-order spectrum');
subplot(2,2,2);
cdfplot(spectrum_plot(index_complete));
set(gca,'yscale','log');
set(gca,'xlim',[0 1]);
xlabel('second-order spectrum');

%pearson
subplot(2,2,3);
hist(pearson_plot(index_complete));
title('Pearson correlation');
set(gca,'xlim',[0 1]);
subplot(2,2,4);
cdfplot(pearson_plot(index_complete));
set(gca,'yscale','log');
xlabel('Pearson correlation');
set(gca,'xlim',[0 1]);

%% check: specturm vs Pearson correlation, perfectly correlated
figure(101);
plot(spectrum_plot(index_complete),pearson_plot(index_complete),'o');
xlabel('Second-order spectrum CDF');
ylabel('Pearson correlation(second-order expansion)');
set(gca,'xlim',[0.4 1],'ylim',[0.4 1]);

%% identify the subgraph with lowest second-order spectrum (higher-order epistasis)
[pearson_sort index_sort]=sort(pearson_plot(index_complete),'ascend');

%candidate list: lowest pearson at second-order expansion
genotype_candidate={};
pearson_candidate=[];
spectrum_candidate=[];
% candidate list: lowest 100 
candidate=100;
for i=1:candidate
pearson_candidate=[pearson_candidate; pearson_plot(index_complete(index_sort(i)))];
spectrum_candidate=[spectrum_candidate; spectrum_plot(index_complete(index_sort(i)))];
genotype_candidate=[genotype_candidate; graph_genotype{index_complete(index_sort(i))}(16,:)];
end

%save
T_candidate= table(genotype_candidate, pearson_candidate, spectrum_candidate);
% writetable(T_candidate,'subgraph_candidate.dat');


%% check: WNWY subgraph (Nick)
% index_WNWY =find(ismember(mutfitness.mut,'WNWY'));
% index_WNWY_HD4=find(index_HD{4}==index_WNWY);
% temp=cell2mat(graph_fitness(index_WNWY_HD4));
% % dlmwrite('fitness_WNWY.txt',temp,'delimiter','\t');
% 
% %plot expansions
% figure(100);
% for i=1:4
% % subplot(2,2,i);
% plot(temp(:,i),temp(:,4),'o');
% hold on;
% end
% x=linspace(min(temp(:,4)),max(temp(:,4)),100);
% y=x;
% plot(x,y,'k-');

%% plot: pick the subgraph with the worst and best second-order expansion (supplementary figure)
figure(101);

%best
subplot(1,2,1);
graph_best=graph_genotype{index_complete(index_sort(end-1))}(16,:);
temp=cell2mat(graph_fitness(index_complete(index_sort(end-1))));
for i=1:4
% subplot(2,2,i);
plot(temp(:,i),temp(:,4),'o');
hold on;
end
% x=linspace(min(temp(:,4)),max(temp(:,4)),100);
% y=x;
x=linspace(-6,2,100);
y=x;
plot(x,y,'k-');
title(graph_best);
legend('1^{st} order','2^{nd} order','3^{rd} order','4^{th} order');
set(gca,'xlim',[-6 2],'ylim',[-6 2],'fontsize',15);
xlabel('Reconstructed fitness');
ylabel('Measured fitness');

%worst
subplot(1,2,2);
graph_worst=graph_genotype{index_complete(index_sort(1))}(16,:);
temp=cell2mat(graph_fitness(index_complete(index_sort(1))));
for i=1:4
% subplot(2,2,i);
plot(temp(:,i),temp(:,4),'o');
hold on;
end
plot(x,y,'k-');
title(graph_worst);
legend('1^{st} order','2^{nd} order','3^{rd} order','4^{th} order');
set(gca,'xlim',[-6 2],'ylim',[-6 2],'fontsize',15);
xlabel('Reconstructed fitness');
ylabel('Measured fitness');

%% plot the correlation for these two subgraphs (supplementary figure)
L=4;
%best
plot(graph_pearson{index_complete(index_sort(end-1))},'-','color','b');
hold on;
%worst
plot(graph_pearson{index_complete(index_sort(1))},'-','color','r');
set(gca,'xlim',[1 L],'ylim',[0 1.1],'fontsize',15);
xlabel('Fourier expansion (order)');
ylabel('Pearson correlation');
legend('CSPA','WLLH');