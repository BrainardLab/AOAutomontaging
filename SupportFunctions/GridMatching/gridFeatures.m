function [f,d,Loc_Index,CNNPos] = gridFeatures(im,useBPFilter,useCNN,CNNParams,w,g,Loc_Index,LocLossPct)

if (nargin < 3)
    useCNN =0;
    CNNParams = [];
end
if (nargin < 2)
    useBPFilter=0;
end

if (~exist('LocLossPct','var'))
    LocLossPct = 0;
end


if (~exist('Loc_Index','var') || isempty(Loc_Index)) %find location if not available
if(useCNN)
    CNNParams.ProbMap.batchsize = 500;
    % Get half patch size
    HalfPatchSize = ceil((CNNParams.PatchSize-1)./2);
    % load in the Net
    load(CNNParams.ProbMap.NetworkPath);
    load(CNNParams.Results.OptimizationPath);

    net = vl_simplenn_move(net, 'cpu');
    net.layers{end}.type = 'softmax';
    ProbParam.PMsigma = OptParam.MaxSigma;
    ProbParam.PMthresh = OptParam.MaxPMthresh;
    ProbParam.ExtMaxH = OptParam.MaxExtMaxH;
    
    
    [CNNPos]= GetConePosSingle(CNNParams,im,net,ProbParam);%find cone locations using CNN
    Loc_Index = sub2ind(size(im), round(CNNPos(:,2)), round(CNNPos(:,1)));
else
    %highpass filter
    if(useBPFilter)
        im = butterworthbpf(im,25,200,4);
    end
    
    %find using local minimum approach
    [im1_mask, X1, Y1] =  findlocalmin(im);
    
    CNNPos =[X1(1:end) Y1(1:end)];
    Loc_Index = sub2ind(size(im), Y1, X1);
    Loc_intensity = im(Loc_Index);
    
    Loc_Index_IntFiltered =  Loc_Index(Loc_intensity > prctile(Loc_intensity,50));%> mean(im1_int)-std(double(im1_int)));
    Loc_Mask_Intfiltered = zeros(size(im));
    Loc_Mask_Intfiltered(Loc_Index_IntFiltered) = 1;
    
    cellsizes1 = findcellsizes(Loc_Index_IntFiltered,Loc_Mask_Intfiltered);
    Loc_Index_SizeFiltered =  Loc_Index_IntFiltered(cellsizes1 > 6);%prctile(cellsizes1,5));%> mean(im1_int)-std(double(im1_int)));
    Loc_Index = Loc_Index_SizeFiltered;
end
else
    CNNPos =[];
end

LN = round(length(Loc_Index)*LocLossPct);
p = randperm(length(Loc_Index));
Loc_Index(p(1:LN)) = [];

Loc_Mask = zeros(size(im));
Loc_Mask(Loc_Index) = 1;


[f, d]=findgridFeature(Loc_Index,Loc_Mask,[],[],[],w,g);
%[f, d]=findgridFeature(I1_filteredsize,im1_filteredsize,im,im1_filtered,im1_mask);
