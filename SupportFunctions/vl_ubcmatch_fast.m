function [matches, scores] = vl_ubcmatch_fast(d1,d2)
%Faster implementation of the vl_ubcmatch function from vl_feat

index1 = 1:size(d1,2);
thresh = 1.5;%default threshold for removing ambiguous matches 

d1 = single(d1);
d2 = single(d2);

%check if any of them are empty
if(size(d1,2) < 1 || size(d2,2) < 1) 
    matches = [];
    scores = [];
    return     
end
MAX_MAT_SIZE = 25000^2;

N1 = size(d1,2);
N2 = size(d2,2);

%calculate all pair-wise distances
% if dotflag

    best = nan(1, N1);
    secondbest = nan(1, N1);
    index2 = nan(1, N1);
    
    chunksize = floor(MAX_MAT_SIZE./N2);

    if chunksize <= N1
        d1chunkinds = uint32(1:chunksize-1:N1);
        if d1chunkinds(end) ~= N1
           d1chunkinds = [d1chunkinds, N1]; 
        end
        tic;
        
        
        for i=1:length(d1chunkinds)-1

            d1chunk = d1(:,d1chunkinds(i):d1chunkinds(i+1));
            
            dist_all = bsxfun(@plus,dot(d1chunk,d1chunk,1)',dot(d2,d2,1)) -2*(d1chunk'*d2);
            
            %find best matches
            [thisbest,thisindex2] = min(dist_all,[],2);

            %remove best matches to find second best matches
            dist_all(sub2ind(size(dist_all),(1:length(thisindex2))',thisindex2))= inf;
            [thissecondbest,I2] = min(dist_all,[],2);

            best(d1chunkinds(i):d1chunkinds(i+1)) = thisbest;
            secondbest(d1chunkinds(i):d1chunkinds(i+1)) = thissecondbest;
            index2(d1chunkinds(i):d1chunkinds(i+1)) = thisindex2;
        end
         toc;

    else
        dist_all = bsxfun(@plus,dot(d1,d1,1)',dot(d2,d2,1))-2*(d1'*d2);

        %find best matches
        [best,index2] = min(dist_all,[],2);

        %remove best matches to find second best matches
        dist_all(sub2ind(size(dist_all),(1:length(index2))',index2))= inf;
        [secondbest,I2] = min(dist_all,[],2);
        index2 = index2';
    end
            
    %is valid if best match is significantly better than second best match
    valid = (best*thresh < secondbest);
    matches = [index1(valid); index2(valid)];
% Alternate approach- memory safe, but ridiculously slow with large
% datasets.
%     [dist_all, index] = pdist2(d2',d1','squaredeuclidean','Smallest',2);
% 
%     best=dist_all(1,:);
%     secondbest=dist_all(2,:);
%     
%     %is valid if best match is significantly better than second best match
%     valid = (best*thresh < secondbest);
%     matches = [index(valid,1); index(valid,2)];
% end


scores = best(valid)';

end


