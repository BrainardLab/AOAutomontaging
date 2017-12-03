function [gg XX YY]= findlocalmin(imageA)
%This code is mostly borrowed from the function 'conecountingfunc.m'
%From the package - 'Cone_Counting_Distribution_Edition_v0p21' written by Robert F Cooper
%Edits by Min Chen

% First pass to find cones and set filter
CutOffinit = 0.6;

Thrshld = 0;  % modification : HT
%verify input validity
ImgDim = size(imageA);


%Begin algorithm
imageA = double(imageA);
fc = imcomplement(imageA);
[M, N] = size(fc);


%FIR filter design
[f1, f2] = freqspace(15, 'meshgrid');
H = ones(15);
fr = sqrt(f1.^2 + f2.^2);
H(fr > CutOffinit) = 0;

window = fspecial('disk', 7);
% window = padarray(window./max(window(:)), [1+(496/2) 1+(496/2)] );
% window = fspecial('disk', 7);
window = window./max(window(:));
h = fwind2(H, window);
fc = imfilter(fc, h, 'replicate', 'same');
% fc = imfilter(fc, h, 0, 'same');

%Morphological markers generation
LocalMins =  imregionalmin(fc, 4);
se = strel('disk', 1, 0);

ConeMark = imdilate(LocalMins, se);

[L, numMark] = bwlabel(ConeMark);
stats = regionprops(L, 'centroid');
X = zeros(numMark, 1);
Y = X;
g = zeros(M, N);

for ii = 1:numMark
    loc = stats(ii).Centroid; %(x, y)
    loc = round(loc); %integral output
    if imageA(loc(2), loc(1)) > Thrshld
        g(loc(2), loc(1)) = 1;
    end
end

g = im2bw(g);
[Y, X] = find(g == 1);

S = [X Y];

% Quicker way to find N-N distance... RFC 06-20-2012
dist_between_pts=squareform(pdist(S)); % Measure the distance from each set of points to the other
max_ident=eye(length(dist_between_pts)).*max(dist_between_pts(:)); % Make diagonal not the minimum for any observation

[minval]=min(dist_between_pts+max_ident,[],2); % Find the minimum distance from one set of obs to another

nmmicronpix  = mean(minval); % Removed the code

conefreqpix = (1./nmmicronpix);
normpowerpix =.5;
CutOffnew = (conefreqpix.*1.2)/normpowerpix;

%Begin algorithm - second time through
ffc = imcomplement(imageA);
[MM, NN] = size(fc);

%FIR filter design - don't need to repeat setup steps as they are the exact
%same. Should save some exec time.
HH = ones(512, 512);
HH(fr > CutOffnew) = 0;
hh = fwind2(HH, window);
ffc = imfilter(ffc, hh, 'replicate', 'same');

%Morphological markers generation
LocalMinsfin = imregionalmin(ffc, 4);
ConeMarkfin = imdilate(LocalMinsfin, se);

[LL, numMarkfin] = bwlabel(ConeMarkfin);
statsfin = regionprops(LL, 'centroid');
XX = zeros(numMarkfin, 1);
YY = XX;
gg = zeros(MM, NN);

for jj = 1:numMarkfin
    loc = statsfin(jj).Centroid; %(x, y)
    loc = round(loc); %integral output
    if imageA(loc(2), loc(1)) > Thrshld
        gg(loc(2), loc(1)) = 1;
    end
end

gg = im2bw(gg);
[YY, XX] = find(gg == 1);

%SS = [XX YY];

% Quicker way to find N-N distance... RFC 06-20-2012
%dist_between_pts=squareform(pdist(SS)); % Measure the distance from each set of points to the other
%max_ident=eye(length(dist_between_pts)).*max(dist_between_pts(:)); % Make diagonal not the minimum for any observation

%[minval, minind]=min(dist_between_pts+max_ident,[],2); % Find the minimum distance from one set of obs to another

%nnmicronfinal = mean(minval.*micronsperpixel);


% Clip edge cones to reduce artifacting
%clipped_coords=coordclip_npoly([YY XX],[2 max(YY)-1],[2 max(XX)-1]);

%XX = clipped_coords(:,2);
%YY = clipped_coords(:,1);

end

