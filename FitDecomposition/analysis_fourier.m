%% Fourier decomposition of complete subgraphs, graph size is 2^4
%last update: 04/06/2015
%codes are adapted from explore_fourier.m
%downstream analysis: analysis_spectrum
clear;
clc;
folder='./analysis_040615/';
%complete graphs
load([folder 'mutfitness']);
load([folder 'graph_index']);
load([folder 'graph_genotype']);
load([folder 'graph_other']);

%saved analysis
% load([folder 'graph_pearson']);
% load([folder 'graph_spectrum']);

%%
A=2;
L=4;
%genotypes of subspace
sequence_total=A^L;
%string format
genotype_str=dec2base(0:sequence_total-1,A); %sorted
%number format
genotype_num=str2num_sequence(genotype_str);
genotype_spin=1-2*genotype_num; 

%% 
%basis of fourier expansion
for i=1:L
    index{i}=combnk(1:L,i);
    coef_num(i)=size(index{i},1);
    q{i}=zeros(coef_num(i),L);
    for j=1:coef_num(i)
        q{i}(j,index{i}(j,:))=1;
        %         q{i}(j,L+1-index{i}(j,:))=1;
    end
end

base=-1;
order_limit=L;

%% analyze Fourier coefficients
total=length(index_HD{4});
graph_spectrum=cell(total,1);
graph_pearson=cell(total,1);
graph_fitness=cell(total,1);

tic;
for n=1:total
    if complete(n)==1 %only analyze complete graphs
        fitness=log(mutfitness.I20fit(graph_index{n}));
        
        %Hadamard-Walsh transform
        y=fwht(fitness,sequence_total,'hadamard');
        % plot(abs(y)/abs(y(1)));
        
        for i=1:order_limit %too costly to compute all orders when L is large (e.g. L=20)
            fq{i}=fouriercoef(base, q{i},fitness, genotype_num);
        end
        
        f0=mean(fitness);
        fourier_all=[];
        fourier_all(1)=1;
        for i=1:order_limit
            fq_norm{i}=fq{i}/f0;
            fourier_all=[fourier_all fq_norm{i}];
            meansquare(i)=mean(fq_norm{i}.^2);
            B(i)=sumsqr(fq_norm{i});
            %     plot(i,fq_norm{i}.^2,'b+');
            %     hold on;
            %     plot(i,meansquare(i),'ro');
            
            % %      plot(i,sumsquare(i),'ro');
            %
            % %     plot(i,abs(fq_norm{i}),'b+');
            % %     hold on;
            % %     meanabs(i)=mean(abs(fq_norm{i}));
            % %     plot(i,meanabs(i),'ro');
        end
        % set(gca,'xlim',[0 order_limit],'yscale','log');
        % xlabel('Order');
        % ylabel('Spectrum');
        
        %contribution from each order
        Fn=B/sum(B);
        for i=1:order_limit
            cdf_Fn(i)=sum(Fn(1:i));
        end
        % plot(cdf_Fn,'-o');
        % set(gca,'xlim',[1 order_limit]);
        % xlabel('Order');
        % ylabel('Spectrum CDF');
        
        % reconstruct fitness using second-order expansion
        for i=1:L
            order=i;
            fitness_fourier(i,:)=reconstructfitness(order,q,fq_norm,genotype_num);
            % plot(fitness,fitness_fourier*f0,'o');
            rankcorr(i)=corr(fitness,fitness_fourier(i,:)'*f0,'type','spearman');
            pearsoncorr(i)=corr(fitness,fitness_fourier(i,:)'*f0,'type','pearson');
            % xlabel('True fitness');
            % ylabel(strcat(num2str(order_limit),'-order Fourier expansion'));
            % ylabel('Second-order Fourier expansion');
            % title(strcat('Spearman=',num2str(rankcorr(i)),',Pearson=',num2str(pearsoncorr(i))));
        end
        
        %save analysis
        graph_spectrum{n}=cdf_Fn;
        graph_pearson{n}=pearsoncorr;
        graph_fitness{n}=fitness_fourier'*f0;
    end
    
end
toc;

%Elapsed time is 451.503689 seconds

%% save
% save('graph_spectrum.mat','graph_spectrum');
% save('graph_pearson.mat','graph_pearson');
% save('graph_fitness.mat','graph_fitness');

%% output data: send to Nick
index_complete=find(complete==1);

genotype_all={};
pearson_all={};
spectrum_all={};
total_list=length(index_complete);
tic;
for i=1:total_list
    genotype_all=[genotype_all; graph_genotype{index_complete(i)}(16,:)]; %16: genotype of the four mutant
    pearson_all=[pearson_all; graph_pearson{index_complete(i)}];
    spectrum_all=[spectrum_all; graph_spectrum{index_complete(i)}];
end
toc;

T_all= table(genotype_all,pearson_all, spectrum_all);
writetable(T_all,'subgraph.dat');
save('subgraph.mat','genotype_all','pearson_all','spectrum_all');

%% save fitness: fourier expansion at all orders

fitness_all={};
total_list=length(index_complete);
tic;
for i=1:total_list
fitness_all=[fitness_all; reshape(graph_fitness{index_complete(i)},[64,1])];
end
toc;

%%
T_fitness= table(genotype_all,fitness_all);
writetable(T_fitness,'subgraph_fitness.dat');



