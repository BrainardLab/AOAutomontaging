The following document contains instructions for reproducing the figures and tables presented in the paper:

"Multi-modal automatic montaging of adaptive optics retinal images"

<Insert Citation on Release>

The software used to produce the data and images can be downloaded at:

https://github.com/BrainardLab/AOAutomontaging

Version 1.0.0 of the software was used to create the figures in the paper.

The data used in the paper can be downloaded at:

Multi-Modal AO Dataset:

https://figshare.com/s/ecaafb98400b0282eaff

Manual Montages:

https://figshare.com/s/2cc83c3cffa0cbbdfd33

Diminishing Overlap Analysis Dataset:

https://figshare.com/s/b548a6989c840fabd832

MATLAB version R2015b was used for processing the data and evaluating the results.

-----Summary of Folder Structure-----

The following is the folder structure you will have by the end of these experiments. Each number is the identifier for the folder, and each period represents a subfolder of the previous number. I.e. (1.1) is a sub folder of (1), while (1.1.2) is a subfolder of (1.1)

(1)OverallProjectDir 

(1.1)AllInputDataDir 

(1.1.1)AOMontagingDataSet (Unzipped from download)

(1.1.2)AOManualMontages (Unzipped from download)

(1.1.3)OverlapAnalysisPairs (Unzipped from download)

(1.2)outputMontageDir (Contents generated using "BatchProcessAOMontages.m")

(1.3)outputAnalysisDir (Contents generated using "CalcAndPlotOverlapSimilarityMetrics.m" and "DiminishingOverlapAnalysis.m")

-----Data Processing Instructions-----

1.) Download the software and data listed above.

2.) Download and install the ToolboxToolbox from:

https://github.com/ToolboxHub/

3.)Copy the file ./Configuration/AOAutomontagingLocalHookTemplate.m into your ToolboxToolbox localHookFolder directory (by default,
% /<Home>/<Matlab>/localHookFolder)

4.)Edit the template file by changing the 6 system paths to match the location of the input data and output locations. (More detailed instructions are provided in the template file).

5.)Rename the file to AOAutomontaging.m, by removing "LocalHookTemplate" from the name.

6.) Call 'tbUse({'AOAutomontaging'})' in Matlab. Make sure there is a notification that the template file above is actually being loaded. If not then please revisit steps 2-5.


7.)Run the script: 

./UtilityAndEvaluationScripts/BatchProcessAOMontages.m

This will automatically montage every available dataset with the transformation and modality settings used in the paper.

*Note this has a high time requirement for completion

-----Figure Reproduction Instructions-----

--Figure 1--

Shown is the temporal arm from the dataset NC_11002_20151006_OS_Images-DONE with and without manual montaging. 

The un-montaged grid (a) was constructed using the nominal location file 11002-20151006.pdf contained with the dataset.

The montaged image (b) was manually constructed and can be found in the Photoshop file NC_11002_20151006_montage_50min_2pc.psd in the dataset.


--Figure 2--

Shown are the following 6 images from the dataset NC_11002_20151006_OS_Images-DONE:
 
NC_11002_123017_OS_confocal_0037_ref_25_lps_8_lbss_8_sr_n_50_cropped_5.tif

NC_11002_123017_OS_confocal_0038_ref_131_lps_8_lbss_8_sr_n_50_cropped_5.tif

NC_11002_123017_OS_split_det_0037_ref_25_lps_8_lbss_8_sr_n_50_cropped_5.tif

NC_11002_123017_OS_split_det_0038_ref_131_lps_8_lbss_8_sr_n_50_cropped_5.tif

NC_11002_123017_OS_avg_0037_ref_25_lps_8_lbss_8_sr_n_50_cropped_5.tif

NC_11002_123017_OS_avg_0038_ref_131_lps_8_lbss_8_sr_n_50_cropped_5.tif

--Figures 3 and 4--

Run the script: 

./UtilityAndEvaluationScripts/ExampleSIFTandRANSACmatch.m

Figure 3 is constructed with the images:

<outputPath>/SIFTmatches_m1.bmp

<outputPath>/SIFTmatches_m2.bmp

<outputPath>/SIFTmatches_m3.bmp


Figure 4 shows the images:

<outputPath>/RANSACMatches_m1.bmp

<outputPath>/RANSACMatches_m2.bmp

<outputPath>/RANSACMatches_m3.bmp

<outputPath>/mosaic_m1.bmp

<outputPath>/mosaic_m2.bmp

<outputPath>/mosaic_m3.bmp

--Figure 5--

If the "Data Processing" section is complete, the automated montage shown in Figure 5 can be found at:

<outputPath>/AutoMontageResults/NC_11048_20151120_OS_Images-DONE/RotationTranslation/AllThree/ref_1_combined_m1.tif

The manual montage shown in the figure can be found at:

<inputPath>/Data/NC_11048_20151120_OS_Images-DONE/NC_11048_20151120_OS_35min_1pc.psd

--Figure 6--

If the "Data Processing" section is complete, run script:

/UtilityAndEvaluationScripts/CalcAndPlotOverlapSimilarityMetrics.m 

The automated montage shown in Figure 6 can be found at:

<outputPath>/CCPlotauto_f10_m1.bmp

The manual montage shown in the figure can be found at:

<outputPath>/CCPlotmanual_f10_m1.bmp

*note the "f10" in the filename may be shifted to "f11" depending on your file system.  

--Figure 7--

Run script:
/UtilityAndEvaluationScripts/DiminishingOverlapAnalysis.m

and then

/UtilityAndEvaluationScripts/PlotDiminishingOverlapAnalysis.m

The three plots shown in Figure 7 will be displayed as MATLAB figures you can save.


-----Table Reproduction Instructions-----

--Table 1 and 3--

Follow the instructions to generate Figure 6. After the script is running, the NCC and NMI numbers in the two tables can be found saved to the files:

<outputPath>/OverlapSimilarityResults.mat

with the variables:

avgCorr_all(f,m,k)  for average of the NCC over all overlapping regions in the dataset

stdDevCorr_all(f,m,k) for standard deviation of the NMI over all overlapping regions in the dataset

avgNMI_all(f,m,k) for standard deviation of the NCC over all overlapping regions in the dataset

stdDevNMI_all(f,m,k) for standard deviation of the NMI over all overlapping regions in the dataset


where 

f ranges from 3 to 13 and indexes the 11 dataset (possibly ranges from 4 to 14 depending on your file system)

m ranges from 1 to 3 and indexes the modality evaluated

m=1 confocal

m=2 split detection

m=2 dark field

k ranges from 1 to 6 and indexes the montaging approach

k=1 is using the manual montage

k=2 is using the automated montage with all three modalities using the rigid transformation model

k=3 is using the automated montage with just the confocal modality using the rigid transformation model

k=4 is using the automated montage with just the split detection modality using the rigid transformation model

k=5 is using the automated montage with just the dark field modality using the rigid transformation model

k=6 is using the automated montage with all three modalities using the translation only transformation model
        

--Table 2--

The discontinuity analysis was performed by manually counting the discontinuities in the manual (Photoshop files) and automated montages (results from the "Data processing" section). And then checking if each discontinuity has as a valid match in the other montage.

(C) Min Chen