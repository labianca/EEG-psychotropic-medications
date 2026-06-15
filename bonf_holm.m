%bonferoni-holm multiple comparisons correction
% from https://www.mathworks.com/matlabcentral/fileexchange/28303-bonferroni-holm-correction-for-multiple-comparisons

function [corrected_p, h]=bonf_holm(pvalues,alpha)
if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvalues<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvalues>1,1)),
        fprintf('WARNING: Some uncorrected p-values are greater than 1.\n');
    end
end
if nargin<2,
    alpha=.05;
elseif alpha<=0,
    error('Alpha must be greater than 0.');
elseif alpha>=1,
    error('Alpha must be less than 1.');
end
s=size(pvalues);
if isvector(pvalues),
    if size(pvalues,1)>1,
       pvalues=pvalues'; 
    end
    [sorted_p sort_ids]=sort(pvalues);    
else
    [sorted_p sort_ids]=sort(reshape(pvalues,1,prod(s)));
end
[dummy, unsort_ids]=sort(sort_ids); %indices to return sorted_p to pvalues order
m=length(sorted_p); %number of tests
mult_fac=m:-1:1;
cor_p_sorted=sorted_p.*mult_fac;
cor_p_sorted(2:m)=max([cor_p_sorted(1:m-1); cor_p_sorted(2:m)]); %Bonferroni-Holm adjusted p-value
corrected_p=cor_p_sorted(unsort_ids);
corrected_p=reshape(corrected_p,s);
h=corrected_p<alpha;
