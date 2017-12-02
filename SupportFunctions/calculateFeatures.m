function [f_all, d_all, h] = calculateFeatures(imageFilename, parallelFlag, pixelScale, featureType, MN, N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%feature parameters
SiftLevel = 55; %The number of levels to use in SIFT, default
ROICropPct = 0; %Sets a percentage crop on the boundaries of the image, where SIFT features are

%stores features for each image
f_all = cell(MN,N);
d_all = cell(MN,N);


h = waitbar(0,'Calculating Sift Features (0%)');
for n=1:N
    if(parallelFlag)
        pixelScale_n = pixelScale(n);
        parfor m = 1:MN
            im = imresize( im2single(imread(char(imageFilename{m,n})) ), pixelScale_n,'bilinear');
            
            if(featureType == 0)
                [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
            elseif(featureType == 1)
                [f1,d1] = gridFeatures(im(:,:,1));
            else
                [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
            end
            
            [f1_crop, d1_crop] = filterSiftFeaturesByROI(im, f1, d1, ROICropPct);
            f_all{m,n} = f1_crop;
            d_all{m,n} = d1_crop;
        end
    else
        for m = 1:MN
            im = imresize( im2single(imread(char(imageFilename{m,n})) ), pixelScale(n),'bilinear');
            if(featureType == 0)
                [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
            elseif(featureType == 1)
                [f1,d1] = gridFeatures(im(:,:,1));
            else
                [f1,d1] = vl_sift(im(:,:,1),'Levels',SiftLevel);
            end
            [f1_crop, d1_crop] = filterSiftFeaturesByROI(im, f1, d1, ROICropPct);
            f_all{m,n} = f1_crop;
            d_all{m,n} = d1_crop;
        end
    end
    waitbar(n/(N),h,strcat('Calculating Sift Features (',num2str(100*n/N,3),'%)'));
end


end

