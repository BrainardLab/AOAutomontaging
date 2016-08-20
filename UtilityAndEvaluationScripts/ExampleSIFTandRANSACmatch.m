%Example script for using SIFT and RANSAC to match a pair of nominal locations with 3 modalities 
%Written by Min Chen (minchen1@upenn.edu)

inputDataDir = 'C:\Users\dontm\Documents\Research\AdaptiveOpticsMosaic\PaperValidationExperiments\BOE_2016\Data\';
outputDir = 'C:\Users\dontm\Downloads\temp\New folder (3)';
datasetDir = 'NC_11002_20151006_OS_Images-DONE\De-identified Data';
MN = 3;%using 3 modalities in this example
N = 2;%two images to match

im1_f = cell(MN,1);
im2_f = cell(MN,1);

im1_f{1}='NC_11002_123017_OS_confocal_0037_ref_25_lps_8_lbss_8_sr_n_50_cropped_5.tif';
im2_f{1}='NC_11002_123017_OS_confocal_0038_ref_131_lps_8_lbss_8_sr_n_50_cropped_5.tif';
im1_f{2}='NC_11002_123017_OS_avg_0037_ref_25_lps_8_lbss_8_sr_n_50_cropped_5.tif';
im2_f{2}='NC_11002_123017_OS_avg_0038_ref_131_lps_8_lbss_8_sr_n_50_cropped_5.tif';
im1_f{3}='NC_11002_123017_OS_split_det_0037_ref_25_lps_8_lbss_8_sr_n_50_cropped_5.tif';
im2_f{3}='NC_11002_123017_OS_split_det_0038_ref_131_lps_8_lbss_8_sr_n_50_cropped_5.tif';

%save SIFT features
f_all = cell(MN,N);
d_all = cell(MN,N);

%read images
im1 = cell(MN,1);
im2 = cell(MN,1);

for m=1:MN
    im1_path=fullfile(inputDataDir,datasetDir,im1_f{m});
    im2_path=fullfile(inputDataDir,datasetDir,im2_f{m});
    im1{m} = imread(im1_path);
    im2{m} = imread(im2_path);
end
%find SIFT features
for m=1:MN
[f1,d1] = vl_sift(im2single(im1{m}),'Levels',55);
f_all{m,1} = f1;
d_all{m,1} = d1;
               
[f2,d2] = vl_sift(im2single(im2{m}),'Levels',55);
f_all{m,2} = f2;
d_all{m,2} = d2;
end
TransType=1;%Use a Rigid example for this example
saveFlag = 1;
%find mataches and RANSAC
sift_mosaic_fast_MultiModal(im1, im2, outputDir,saveFlag,f_all(:,1),d_all(:,1),f_all(:,2),d_all(:,2),TransType)