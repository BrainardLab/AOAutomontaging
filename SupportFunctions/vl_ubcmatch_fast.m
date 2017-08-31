function [matches, scores] = vl_ubcmatch_fast(d1,d2)
%Faster implementation of the vl_ubcmatch function from vl_feat

index1 = 1:size(d1,2);
thresh = 1.5;%default threshold for removing ambiguous matches 

d1 = single(d1);
d2 = single(d2);

%calculate all pair-wise distances

dist_all = pdist2(d1',d2','squaredeuclidean','Smallest',2);

best=dist_all(:,1);
secondbest=dist_all(:,2);

dist_all = bsxfun(@plus,dot(d1,d1,1)',dot(d2,d2,1))-2*(d1'*d2);

%find best matches
[best,index2] = min(dist_all,[],2);

%remove best matches to find second best matches
dist_all(sub2ind(size(dist_all),(1:length(index2))',index2))= inf;
[secondbest,I2] = min(dist_all,[],2);

%is valid if best match is significantly better than second best match
valid = (best*thresh < secondbest);

%output results
matches = [index1(valid); index2(valid)'];
scores = best(valid)';

end


