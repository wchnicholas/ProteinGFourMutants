function epistasis_type=epistasis_classify(fitness_quad,genotype_quad)
  %classify types of pairwise epistasis (algorithm in Greene and Crona, 2014)
  %class 1: magnitude; class 2: sign; class 3: reciprocal sign
  L=size(genotype_quad,2); %genome length
  
  %fitness, rank in ascending order
  [fitness_ascend, index_ascend]=sort(fitness_quad);
  %index of genotypes,rank in ascending order
  genotype_ascend=genotype_quad(index_ascend,:);
  
  %1) magnitude epistasis: if rank#1 and #4 are not nearest neighbors
  %TBD: separate magnitude epistasis into synergistic and antagonistic, assuming log(fitness) is additive
  if pdist2(genotype_ascend(1,:),genotype_ascend(4,:),'hamming')==2/L
      epistasis_type=1;
      %2) reciprocal sign epistasis: if #3 and #4 are not nearest neighbors
  else if pdist2(genotype_ascend(3,:),genotype_ascend(4,:),'hamming')==2/L
          epistasis_type=3;
          %3) sign epistasis: if #2 and #4 are not nearest neighbors
      else
          epistasis_type=2;
      end
  end
  
end