function [graph_index, complete]= findinter_index(all_genotype, graph_genotype, index_HD, index_WT, Hdistance)
%genotype: a list of genotypes
%Hamming distance to WT specified in Hdistance

total=length(index_HD{4}); %quadruple mutants
graph_index=cell(total,1);
complete=ones(total,1);

for i=1:total
    graph_index{i}=zeros(16,1);
    graph_index{i}(1)=index_WT;
    graph_index{i}(16)=index_HD{4}(i);

    for j=2:15
        distance=Hdistance(j);
        %Note: limit the search to a subset of genotype with the same HD
        index=find(ismember(all_genotype(index_HD{distance}),graph_genotype{i}(j,:))); 
        if isempty(index)
            graph_index{i}(j)=0;
            complete(i)=0;
        else
        graph_index{i}(j) =index_HD{distance}(index);
        end
    end
    
end

end