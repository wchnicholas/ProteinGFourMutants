%% analyze properties of subgraphas
%updated: 11/19/2015
clear;
clc;
%complete subgraphs identified in analysis_040615
folder='./analysis_040615/';
load([folder 'mutfitness']);
load([folder 'graph_index']);
load([folder 'graph_genotype']);
load([folder 'graph_other']);


%% genotypes of subgraph: binary sequence: adapted from constructNK.m
A=2;
L=4;
sequence_total=A^L; %total number of genotypes in the sequence space
% generate a list of genotypes (enumerate the entire sequence space)
%string format
genotype_str=dec2base(0:sequence_total-1,A); %sorted
%number format
genotype_num=str2num_sequence(genotype_str);

%% find subgraphs where quadruple mutants are the single peak
%compare to Nick's: the number should be 29
%1)quadruple mutants is a local peak in the subgraph
%2)exclude subgraphs with more than one peak 
%(i.e. quadruple mutant is the global peak)
%consider: exclude subgraphs with adjacent lethal mutants?

total=length(index_HD{4});
count=0;
subgraph_adapt=[];
tic;
for n=1:total
    if complete(n)==1 %only analyze complete graphs
        fitness_subgraph=log(mutfitness.I20fit(graph_index{n}));
        %find HD=3 mutants
        index_triple=find(Hdistance==3);        
        fitness_triple=fitness_subgraph(index_triple);
        fitness_quadruple=fitness_subgraph(16);
        if fitness_quadruple>max(fitness_triple) %if quadruple mutant is a local peak in the subgraph
            [peaknum,basin_size] = localoptima(fitness_subgraph,genotype_str); %localoptima.m
            if peaknum==1
                count=count+1;
                subgraph_adapt(count)=n;
            end
        end
    end
end
toc;

%% analysis
%1)count accessible pathways: adapted from accessiblepath.m
for i=1:length(subgraph_adapt)
     fitness_subgraph=log(mutfitness.I20fit(graph_index{subgraph_adapt(i)}));

     %find all the paths in this subgraph (permutation, total # of paths=Mdistance!)
     %record the index of intermediate sequences and fitness trajectory along each path
     Mdistance=L;
     path_total=factorial(Mdistance);
     v=1:Mdistance;
     mutation_order=perms(v);
     index_order=zeros(path_total,Mdistance-1);
     fitness_order=zeros(path_total,Mdistance-1);
     %start at WT (0000)
     index_start=1;
     fitness_start=fitness_subgraph(index_start);
     
     %enumerate all possible paths   
     for path=1:path_total %path
         index_order(path,1)=index_start;
         fitness_order(path,1)=fitness_subgraph(index_start);
         for k=1:Mdistance %mutation
             index_order(path,k+1)=findmutant(index_order(path,k),mutation_order(path,k),genotype_num); %findmutant.m
             fitness_order(path,k+1)=fitness_subgraph(index_order(path,k+1));
         end
     end
     %calculate the fraction of monotonically increasing path
     accessible=zeros(path_total,1);
     for n=1:path_total
         if all(diff(fitness_order(n,:))>0)
             accessible(n)=1;
         end
     end
     accessible_sum(i)=sum(accessible);
     accessible_frac(i)=sum(accessible)/path_total;
     
end
%% bar plot
%compare to Nick's: identical
accessible_sum_descend=sort(accessible_sum,'descend');
bar(accessible_sum_descend);

%% analysis
%2)enumerate pairwise epistasis, calculate the fraction of different types 
epistasis_count=zeros(length(subgraph_adapt),3);
for i=1:length(subgraph_adapt)
     fitness_subgraph=log(mutfitness.I20fit(graph_index{subgraph_adapt(i)}));
     count=0;
     epistasis_type=[];
     for j=1:sequence_total %enumerate genotypes
         index_start=j;
         mutation_pair=combnk(1:L,2); %enumerate pairs of single mutations
         for k=1:length(mutation_pair)
             index_single1=findmutant(index_start,mutation_pair(k,1),genotype_num);
             index_single2=findmutant(index_start,mutation_pair(k,2),genotype_num);
             index_double=findmutant(index_start,mutation_pair(k,:),genotype_num);
             index_quad=[index_start,index_single1,index_single2,index_double];
             fitness_quad=fitness_subgraph(index_quad);
             genotype_quad=genotype_num(index_quad,:);
             count=count+1;
             epistasis_type(count)=epistasis_classify(fitness_quad, genotype_quad); %epistasis_classify.m
         end
     end
     for m=1:3
         epistasis_count(i,m)=length(find(epistasis_type==m))/4; %divide by 4, because every quad is counted 4 times. total=16*6/4=24
     end
end

%% ruggedness score
epistasis_total=24;
for i=1:length(subgraph_adapt)
    ruggedness(i)=(epistasis_count(i,2)+epistasis_count(i,3)*2)/epistasis_total;
end

%% scatter: some dots are overlapping
plot(ruggedness,path_total-accessible_sum,'o',...
    'MarkerFaceColor','b');
xlabel('Ruggedness (Epistasis)');
ylabel('Ruggedness (Inaccessible paths)');
set(gca,'fontsize',15','xlim',[0 1],'ylim',[10 24]);
box off;

%% scatter: mark dots with >1 subgraph
%source:http://stackoverflow.com/questions/13778799/scatterplot-visualize-the-same-points-in-matlab
A=[ruggedness' path_total-accessible_sum'];
[Auniq,~,IC] = unique(A,'rows');
cnt = accumarray(IC,1);
scatter(Auniq(:,1), Auniq(:,2), (8+4*(cnt>1)).^2, 'LineWidth',2,'markeredgecolor','k'); % make the ones where we'll put a number inside a bit bigger
for ii=1:numel(cnt)
    if cnt(ii)>1
        text(Auniq(ii,1),Auniq(ii,2),num2str(cnt(ii)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize', 8, 'color','k');
    end
end
xlabel('Ruggedness (pairwise epistasis)');
ylabel('Inaccessible direct paths');
set(gca,'fontsize',15','xlim',[0 0.8],'ylim',[10 24]);
box off;

%% correlation
[corr_ruggedness,pvalue]=corr(ruggedness',path_total-accessible_sum','type','pearson')

% corr_ruggedness =
%     0.6585
% pvalue =
%    1.0295e-04
