I1 = im2gray(imread('../Demos/Multi-Modal/CS_13213_20150302_confocal_OS_0070_ref_12_lps_8_lbss_8_ffr_n_50_cropped_5.tif'));
I2 = im2gray(imread('../Demos/Multi-Modal/CS_13213_20150302_confocal_OS_0071_ref_51_lps_8_lbss_8_ffr_n_50_cropped_5.tif'));
%Find the corners.

points1 = detectSIFTFeatures(I1,NumLayersInOctave=5);
points2 = detectSIFTFeatures(I2,NumLayersInOctave=5);
%Extract the neighborhood features.

[features1,valid_points1] = extractFeatures(I1,points1);
[features2,valid_points2] = extractFeatures(I2,points2);
%Match the features.

indexPairs = matchFeatures(features1,features2,MatchThreshold=100);
%Retrieve the locations of the corresponding points for each image.

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
%Visualize the corresponding points. You can see the effect of translation between the two images despite several erroneous matches.

figure; 
showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);