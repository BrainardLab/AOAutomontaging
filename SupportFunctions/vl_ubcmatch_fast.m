function [matches, scores] = vl_ubcmatch_fast(d1,d2)
%Faster implementation of the vl_ubcmatch function from vl_feat

thresh = 1.5;%default threshold for removing ambiguous matches 

d1_single = single(d1);
d2_single = single(d2);

%calculate all pair-wise distances
dist_all = bsxfun(@plus,dot(d1_single,d1_single,1)',dot(d2_single,d2_single,1))-2*(d1_single'*d2_single);

%find best matches
[best,index2] = min(dist_all,[],2);

%remove best matches to find second best matches
dist_all(sub2ind(size(dist_all),(1:length(index2))',index2))= inf;
[secondbest,I2] = min(dist_all,[],2);

%is valid if best match is significantly better than second best match
valid = (best*thresh < secondbest);

%output results
index1 = 1:size(d1,2);
matches = [index1(valid); index2(valid)'];
scores = best(valid)';

end


