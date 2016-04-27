%analyze the use of bypass to circumvent reciprocal sign epistasis
%updated: 11/19/2015
clear;
clc;

%load original data
%fitness
load analysis_040615/mutfitness
%genotype_int
load bypass/genotype_int
%reciprocal epistasis
load bypass/index_reciprocal_1e5

%% bypass via detour
sample_reciprocal=length(index_reciprocal);
% sample_reciprocal=10;
tic;
for i=1:sample_reciprocal
    %for the reciprocal epistasis, find the starting point(rank#3) and end point (rank#4)
    genotype_sample=double(genotype_int(index_reciprocal{i},:));
    fitness_sample=mutfitness.I20fit(index_reciprocal{i}); %this is not ordered
    %fitness, rank in ascending order
    [fitness_ascend, index_ascend]=sort(fitness_sample);
    %index of genotypes,rank in ascending order
    genotype_ascend=genotype_sample(index_ascend,:);
    %starting point, rank#3 (000)
    genotype_000=genotype_ascend(3,:);
    %end point, rank#4 (011)
    genotype_011=genotype_ascend(4,:);
    %identify the position of two mutations
    mutation_pos=find(genotype_011-genotype_000);
    genotype_001=genotype_ascend(1,:); %here the assignment of 001 or 010 does not matter (symmetric)
    genotype_010=genotype_ascend(2,:);
    %fitness, reordered 000,001,010,011
    index_000=index_reciprocal{i}(index_ascend(3));
    index_011=index_reciprocal{i}(index_ascend(4));
    index_001=index_reciprocal{i}(index_ascend(1));
    index_010=index_reciprocal{i}(index_ascend(2));
    index_quad=[index_000,index_001,index_010,index_011];
    fitness_quad=mutfitness.I20fit(index_quad);
    
    %enumerate all possible mutations on an additional third position
    position_all=1:4;
    position_additional=setdiff(position_all,mutation_pos); %sites for additional mutation
    for j=1:length(position_additional) 
        %identify mutation spectrum of residue
        site=position_additional(j);
        residue=genotype_000(site);
        temp=1:1:20;
        mutationspec{j}=temp(temp~=residue);
    end
    
    %for each third mutation, identify if a detour bypass is possible
    count_all(i)=0; %all possible routes (excluding missing genotypes)
    count_bypass(i)=0;
    for j=1:length(position_additional) %enumerate position (4-2=2 in total)
        for aa=1:length(mutationspec{j}) %enumerate mutation spectrum (amino acid)
            site=position_additional(j);
            %find genotype and fitness of (100,101,110,111)
            genotype_100=genotype_000;
            genotype_100(site)=mutationspec{j}(aa);
            genotype_111=genotype_011;
            genotype_111(site)=mutationspec{j}(aa);
            %intermediates
            genotype_101=genotype_001;
            genotype_101(site)=mutationspec{j}(aa);
            genotype_110=genotype_010;
            genotype_110(site)=mutationspec{j}(aa);
            %index
            index_100=find(ismember(mutfitness.mut,int2aa(genotype_100)));
            index_111=find(ismember(mutfitness.mut,int2aa(genotype_111)));
            index_101=find(ismember(mutfitness.mut,int2aa(genotype_101)));
            index_110=find(ismember(mutfitness.mut,int2aa(genotype_110)));
            %exclude if genotypes are missing
            if ~isempty(index_100) && ~isempty(index_111) && ~isempty(index_101) && ~isempty(index_110)
                count_all(i)=count_all(i)+1;
                %fitness
                index_quad_new=[index_100,index_101,index_110,index_111];
                fitness_quad_new=mutfitness.I20fit(index_quad_new);
                %fitness, rank in ascending order
                [fitness_ascend_new, index_ascend_new]=sort(fitness_quad_new);
                
                %criterion for a successful detour bypass
                %1)100>000 (detour step is beneficial);
                %2)111<011 (loss);
                %3)111 is rank#4 among (100,101,110,111)
                %4)100 is not rank #3 among (100,101,110,111)
                if fitness_quad_new(1)>fitness_quad(1) && fitness_quad_new(4)<fitness_quad(4) && index_ascend_new(end)==4 &&index_ascend_new(end-1)~=1
                    count_bypass(i)=count_bypass(i)+1;
                end
            end
        end
    end
end
toc;

%241 samples (out of 1000 pairwise epistasis)
%Elapsed time is 265.221624 seconds.

%% bypass via conversion

sample_reciprocal=length(index_reciprocal);
% sample_reciprocal=1;
tic;
for i=1:sample_reciprocal
    %for the reciprocal epistasis, find the starting point(rank#3) and end point (rank#4)
    genotype_quad=double(genotype_int(index_reciprocal{i},:));
    fitness_quad=mutfitness.I20fit(index_reciprocal{i});
    %fitness, rank in ascending order
    [fitness_ascend, index_ascend]=sort(fitness_quad);
    %index of genotypes,rank in ascending order
    genotype_ascend=genotype_quad(index_ascend,:);
    %starting point, rank#3 (00)
    genotype_00=genotype_ascend(3,:);
    fitness_00=fitness_ascend(3);
    %end point, rank#4 (11)
    genotype_11=genotype_ascend(4,:);
    fitness_11=fitness_ascend(4);
    %identify the position of two mutations
    mutation_pos=find(genotype_11-genotype_00);

    %enumerate all possible mutations at these 2 positions
    for j=1:length(mutation_pos) 
        %identify mutation spectrum of residue
        site=mutation_pos(j);
        %exclude 01, 10
        temp=1:1:20;
        mutationspec{j}=temp(temp~=residue);
    end
    
    %for each mutation, identify if a conversion bypass is possible
    count_all_conversion(i)=0; %all possible routes (excluding missing genotypes)
    count_bypass_conversion(i)=0;
    for j=1:length(mutation_pos) %enumerate position (2 in total)
        for aa=1:length(mutationspec{j}) %enumerate mutation spectrum (amino acid)
            %the position with conversion
            site=mutation_pos(j);
            %the other position
            site_other=mutation_pos(mutation_pos~=mutation_pos(j));
            
            %find genotype and fitness of (20,21)
            genotype_20=genotype_00;
            genotype_20(site)=mutationspec{j}(aa);
            genotype_21=genotype_20;
            residue_other=genotype_11(site_other);
            genotype_21(site_other)=residue_other;

            %index
            index_20=find(ismember(mutfitness.mut,int2aa(genotype_20)));
            index_21=find(ismember(mutfitness.mut,int2aa(genotype_21)));
            fitness_20=mutfitness.I20fit(index_20);
            fitness_21=mutfitness.I20fit(index_21);
    
            %exclude if genotypes are missing
            if ~isempty(index_20) && ~isempty(index_21) 
                count_all_conversion(i)=count_all_conversion(i)+1;
                %criterion for a successful conversion bypass
                %1)20>00 (conversion step is beneficial)
                %2)21>20
                %3)11>21
                if fitness_20>fitness_00 && fitness_21>fitness_20 && fitness_11>fitness_21
                    count_bypass_conversion(i)=count_bypass_conversion(i)+1;
                end
            end
        end

    end
end
toc;

%1e3
%Elapsed time is 135.637023 seconds.

%% time: 1e5 pairwise (~2e4 reciprocal sign samples)
%Elapsed time is 22411.735763 seconds.
%Elapsed time is 11399.727887 seconds.

%% save
% save('./bypass/bypassfreq_1e5','count_all','count_bypass','count_all_conversion','count_bypass_conversion');

%% plot: see analysis_bypass_plot

% %% calculate the frequency of successful bypass among all possible trials
% %detour 
% %if no missing genotype, count_all=19*(4-2)=38
% freq_bypass=count_bypass./count_all;
% %exclude entries if freq_bypass=NaN (because count_all=0)
% freq_bypass=freq_bypass(~isnan(freq_bypass));
% 
% %conversion
% %if no missing genotype, count_all_conversion=19*2=38
% %NOTE: count_all_conversion should be 18*2=36, 12/4/2015
% freq_bypass_conversion=count_bypass_conversion./count_all_conversion;
% %exclude entries if freq_bypass=NaN (because count_all=0)
% freq_bypass_conversion=freq_bypass_conversion(~isnan(freq_bypass_conversion));
% 
% %% fit to poisson distribution
% [lambdahat_detour,lambdaci_detour] = poissfit(freq_bypass)
% [lambdahat_conversion,lambdaci_conversion] = poissfit(freq_bypass_conversion)
% 
% %outputs
% % lambdahat_detour =
% %     0.0070
% % lambdaci_detour =
% %     0.0059
% %     0.0081
% % 
% % lambdahat_conversion =
% %     0.0321
% % lambdaci_conversion =
% %     0.0297
% %     0.0345
% 
% %% compare two mechanism
% freq_max=max(max(freq_bypass),max(freq_bypass_conversion));
% subplot(1,2,1);
% hist(freq_bypass);
% % set(gca,'xlim',[0 freq_max]);
% set(gca,'ylim',[0 2.1e4]);
% xlabel('Bypass availability');
% ylabel('Reciprocal epistasis');
% title('Detour');
% % hold on;
% % x = 0:0.001:freq_max;
% % y = poisspdf(x,lambdahat_detour)*length(freq_bypass);
% % plot(x,y,'--r')
% 
% subplot(1,2,2);
% hist(freq_bypass_conversion);
% % set(gca,'xlim',[0 freq_max]);
% set(gca,'ylim',[0 2.1e4]);
% xlabel('Bypass availability');
% ylabel('Reciprocal epistasis');
% title('Conversion');
% % hold on;
% % x = 0:0.001:freq_max;
% % y = poisspdf(x,lambdahat_conversion)*length(freq_bypass_conversion);
% % plot(x,y,'--r')

