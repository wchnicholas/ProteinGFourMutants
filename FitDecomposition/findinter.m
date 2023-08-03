function [genotype]=findinter(WT,mut)
%find intermediate genotypes
L=4;
total=2^L;
genotype_str=dec2base(0:total-1,2);
genotype_num=str2num_sequence(genotype_str);

% single=nchoosek(L,1);
% double=nchoosek(L,2);
% triple=noochosek(L,3);

for i=1:total
    for j=1:L
    if genotype_num(i,j)==0
        genotype(i,j)=WT(j);
    else
        genotype(i,j)=mut(j);
    end
    end
end
    

end