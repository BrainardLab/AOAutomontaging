function [f_all, d_all, Loc_Index] = calculateFeatures_Debug(imageFilename, parallelFlag, pixelScale, featureType, MN, N,BPFilterFlags,CNNFlags,w,g,Loc_Index,LocLossPct)
%Calculate features for all images listed in imageFilename

%feature parameters
SiftLevel = 55; %The number of levels to use in SIFT, default
ROICropPct = 0; %Sets a percentage crop on the boundaries of the image, where SIFT features are

%stores features for each image
f_all = cell(MN,N);
d_all = cell(MN,N);

if (featureType == 0)
    FeatureName = 'SIFT';
elseif (featureType == 1)
    FeatureName = 'Constellation';
end

if ~exist('BPFilterFlags','var')
  BPFilterFlags=zeros(1,MN);
end

if ~exist('CNNFlags','var')
    CNNFlags=zeros(1,MN);
else
    % choose dataset with already trained cnn and detection parameters
    CNNParams=cell(1,MN); 
    DataSet = 'confocal';
    CNNParams{1} = get_parameters_Cone_CNN(DataSet);
    DataSet = 'split detector';
    CNNParams{2} = get_parameters_Cone_CNN(DataSet);
    CNNParams{3} = [];
    MatConvNetPath = fullfile('C:\Users\dontm\Documents\OpenSourceLibraries\Matlab\matconvnet-1.0-beta25');
    run(fullfile(MatConvNetPath,'matlab','vl_setupnn.m'))
end

if(sum(CNNFlags)>0) %Can't parallel process if using CNN
    parallelFlag = 0;
end

if (~exist('Loc_Index','var') || isempty(Loc_Index))
    Loc_Index = cell(MN,N);
end

%TODO: Right now the direct lcoation input is not modality specific, same
%values gets set for all modalities

%h = waitbar(0,['Calculating ' FeatureName ' Features (0%)']);
for n=1:N
    if(parallelFlag)
        pixelScale_n = pixelScale(n);
        parfor m = 1:MN
            if ~isempty(imageFilename{m,n}) % If this file is blank, then that means we don't have valid information for it- skip it.
                im = imresize( im2single(imread(char(imageFilename{m,n})) ), pixelScale_n,'bilinear');
                
                if(featureType == 0)
                    [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
                elseif(featureType == 1)
                        [f1,d1] = gridFeatures(im(:,:,1),BPFilterFlags(m),0,w,g,Loc_Index{m,n},LocLossPct);%use bandpassfilter if split, can't use CNN if parallel processing
                else
                    [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
                end
                
                [f1_crop, d1_crop] = filterSiftFeaturesByROI(im, f1, d1, ROICropPct);
                f_all{m,n} = f1_crop;
                d_all{m,n} = d1_crop;
            end
        end
    else
        for m = 1:MN
            if ~isempty(imageFilename{m,n}) % If this file is blank, then that means we don't have valid information for it- skip it.
                im = imresize( im2single(imread(char(imageFilename{m,n})) ), pixelScale(n),'bilinear');
                if(featureType == 0)
                    [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
                elseif(featureType == 1)
                    [f1,d1,Loc_Index_m_n] = gridFeatures(im(:,:,1),BPFilterFlags(m),CNNFlags(m),CNNParams{m},w,g,Loc_Index{m,n},LocLossPct);%use CNN
                else
                    [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
                end
                [f1_crop, d1_crop] = filterSiftFeaturesByROI(im, f1, d1, ROICropPct);
                f_all{m,n} = f1_crop;
                d_all{m,n} = d1_crop;
                Loc_Index{m,n} = Loc_Index_m_n;
            end
        end
    end
%    waitbar(n/(N),h,['Calculating ' FeatureName ' Features (' num2str(100*n/N,3) '%)']);
end


end

