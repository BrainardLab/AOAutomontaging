function matches=matchGridFeatures(d1,d2,d1_int,d2_int)

N1 = size(d1,2);
N2 = size(d2,2);
dist_all = zeros(N1,N2);
thresh = 1.05;
thresh2 = .2;
index1 = 1:size(d1,2);

d1=single(d1);
d2=single(d2);

% for i = 1:N1
%     for j =  1:N2
%          dist_all(i,j) = sum(d1(:,i) & d2(:,j))/(sum(d1(:,i))+sum(d2(:,j))); 
%     end
% end

% Linear alegbra-based matching approach is much faster, especially in MATLAB.
sumd1 = sum(d1);
sumd2 = sum(d2);

dist_all = d1' * d2 ./ bsxfun(@plus,sumd1',sumd2);

        [best,index2] = max(dist_all,[],2);

        %dist_all(isnan(dist_all))=0;
        %dist_all(isinf(dist_all))=0;

        
        %find best matches
        %dist_all(sub2ind(size(dist_all),(1:length(index2))',index2))= 0;
        %[secondbest,I2] = max(dist_all,[],2);
        index2 = index2';
        %valid = (best*thresh > secondbest);
        valid = (best > thresh2);
        matches = [index1(valid); index2(valid)];
        
%        savegridfeaturePairs(d1,d2,d1_int,d2_int,matches,8);
        
