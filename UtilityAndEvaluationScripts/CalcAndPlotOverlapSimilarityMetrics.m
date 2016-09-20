%Script for calculating NCC and NMI similarity metrics in overlapping
%regions of completed montages
%Written by Min Chen (minchen1@upenn.edu)

%load manual and automated montage result folders 
manualDirBase = getpref('AOAutomontaging','inputManualDataDir');
autoDirBase = getpref('AOAutomontaging','outputMontageDir');
outDirBase= getpref('AOAutomontaging','outputAnalysisDir');

manualFileFolders = dir(manualDirBase);
autoFileFolders = dir(autoDirBase);

FN = length(manualFileFolders);
MN = 3;
KN = 6;
avgCorr_all = zeros(FN,MN,KN);
avgNMI_all = zeros(FN,MN,KN);

stdDevCorr_all = zeros(FN,MN,KN);
stdDevNMI_all= zeros(FN,MN,KN);

for f = 1:FN%Do for each dataset
    f
    subFileFolder = manualFileFolders(f).name;
    if strcmp(subFileFolder,'.')||strcmp(subFileFolder,'..')
        continue;
    end
    
    %Set folders we pull the montages from
    autoDir = fullfile(autoDirBase,subFileFolder,'RotationTranslation','AllThree');
    autoConfocalDir = fullfile(autoDirBase,subFileFolder,'RotationTranslation','Confocal');
    autoDFDir = fullfile(autoDirBase,subFileFolder,'RotationTranslation','DarkField');
    autoSplitDir = fullfile(autoDirBase,subFileFolder,'RotationTranslation','SplitDetection');
    autoTransOnlyDir=fullfile(autoDirBase,subFileFolder,'TranslationOnly','AllThree');
    manualDir = fullfile(manualDirBase,subFileFolder,'montaged_individual');
    modality = ['confocal' 'split_det' 'avg'];
    %Calculate Metric for each Modality (m) and Each config (k)
    for m = 1:MN
        m
        switch m %set search string depending on the modality
            case 1
                srchstr = '\*confocal*.tif';
            case 2
                srchstr = '\*split_det*.tif';
            case 3
                srchstr = '\*avg*.tif';
        end
        
        for k = 1:KN
            k
            if(k == 3 && m ~= 1)%k=3 is Confocal Only, skip other modalities
                continue;
            end
            if(k == 4 && m ~= 2)%k=4 is Split Only, skip other modalities
                continue;
            end
            if(k == 5 && m ~= 3)%k=5 is DF Only, skip other modalities
                continue;
            end
            
            %Set input directory depending k
            if (k == 1) %manual
                ImageDir = manualDir;
            elseif (k == 2)%Auto using all 3 modalities
                ImageDir = autoDir;
            elseif (k == 3)%Auto Using Confocal Only
                ImageDir = autoConfocalDir;
            elseif (k == 4)%Auto Using Split Only
                ImageDir = autoSplitDir;
            elseif (k == 5)%Auto Using DF Only
                ImageDir = autoDFDir;
            elseif (k == 6)%Auto Using All 3 Modalities, but Trans Only
                ImageDir = autoTransOnlyDir;
            end
            
            %search Modality to evaluate
            ImageAllfiles = dir(strcat(ImageDir,srchstr));
            ImageAllfiles = {ImageAllfiles.name};
            
            N = size(ImageAllfiles,2);
            imSize = size(imread(fullfile(ImageDir, ImageAllfiles{1})));
            CorrelationPlot = zeros(imSize(1:2));
            NMIPlot = zeros(imSize(1:2));
            
            OverlapCount = zeros(imSize(1:2));
            FullMask  = zeros(imSize(1:2));
            
            im_n = cell(1,N);
            im_n_mask = cell(1,N);
            ref = zeros(1,N);
            
            %Check through montage first
            for i = 1:N 
                i
                if(k ~= 1)%check which reference frame
                    
                    fname = ImageAllfiles{i};
                    indexStart=strfind(fname,'aligned_to_ref');
                    %find how many digits
                    suffix = fname(indexStart:end);
                    indexEnd=indexStart+strfind(suffix,'_m')-1;
                    refstr = fname(indexStart+14:indexEnd-1);
                    ref(i) = str2double(refstr);
                end
                
                %Check image type (RGB, Gray scale, etc.)
                if(length(imSize) == 3)
                    if(imSize(3) == 2)
                        im1 = imread(fullfile(ImageDir, ImageAllfiles{i}));
                        im1 = im1(:,:,1);
                    else
                        im1 = rgb2gray(imread(fullfile(ImageDir, ImageAllfiles{i})));
                    end
                elseif(length(imSize) == 2)
                    im1 = imread(fullfile(ImageDir, ImageAllfiles{i}));
                end
                
                im1(im1 == 255) = 0;
                im_n{i} = im1;
                im_n_mask{i} = im1 > 0;
                FullMask(im_n_mask{i}) = 1;
            end
            
            correlationsPairs = nan(N,N);
            mutualInformationPairs = nan(N,N);
            overlapPairs = nan(N,N);
            
            
            %Go through montage pair by pair
            for i = 1:N
                for j = (i+1):N
                    imIntersect = im_n_mask{i} & im_n_mask{j};
                    %make sure there is an intersection of >1 pixel and
                    %both images are in the same reference frame
                    if(sum(imIntersect(:)) > 1 && (ref(i) == ref(j)))
                        
                        corr_ij = corr2(im_n{i}(imIntersect),im_n{j}(imIntersect));
                        if(~isnan(corr_ij))%make sure the pair can actually be compared
                            correlationsPairs(i,j) = corr_ij;
                            CorrelationPlot(imIntersect) = CorrelationPlot(imIntersect) + corr_ij;
                            
                            nmi_ij = nmi(double(im_n{i}(imIntersect)),double(im_n{j}(imIntersect)));
                            mutualInformationPairs(i,j) = nmi_ij;
                            NMIPlot(imIntersect) = NMIPlot(imIntersect) + nmi_ij;
                            
                            
                            OverlapCount = OverlapCount + imIntersect;
                            overlapPairs(i,j) = sum(sum(imIntersect));
                        end
                    end
                end
            end
            
            %cal and save
            avgCorr = sum(CorrelationPlot(:))/sum(OverlapCount(:));
            avgNMI = sum(NMIPlot(:))/sum(OverlapCount(:));
            
            allCorrPairs = correlationsPairs(~isnan(correlationsPairs));
            allOverlapPairs = overlapPairs(~isnan(overlapPairs));
            allNMIPairs = mutualInformationPairs(~isnan(correlationsPairs));
            
            stdDevCorr = sqrt(sum(((allCorrPairs-avgCorr).^2).*allOverlapPairs)/sum(allOverlapPairs));
            stdDevNMI = sqrt(sum(((allNMIPairs-avgNMI).^2).*allOverlapPairs)/sum(allOverlapPairs));
            
            avgCorr_all(f,m,k) = avgCorr;
            avgNMI_all(f,m,k) = avgNMI;
            stdDevCorr_all(f,m,k) = stdDevCorr;
            stdDevNMI_all(f,m,k) = stdDevNMI;

            plotOn = 1;%set this flag if you want to plot metric maps
            if(plotOn)
                
                if(k==1)
                    rater='manual';
                elseif(k==2)
                    rater='auto';
                elseif(k==3)
                    rater='autoConfocal';
                elseif(k==4)
                    rater='autoSplit';
                elseif(k==5)
                    rater='autoDF';
                elseif(k==6)
                    rater='autoTransOnly';
                end
                
                figID1=figure(1); % avg correlation plot
                q = OverlapCount;
                q(q==0) = 1;
                imshow(CorrelationPlot./q - ~FullMask(:,:,1)*.1)
                colormap('jet')
                caxis([-.1 1])
                
                set(gca,'LooseInset',get(gca,'TightInset'));
                saveMatchesName=['CCPlot' rater '_f' num2str(f) '_m' num2str(m) '.bmp'];
                saveas(figID1,fullfile(outDirBase,saveMatchesName))
                
                figID2=figure(2);% avg nmi plot
                q = OverlapCount;
                q(q==0) = 1;
                imshow(NMIPlot./q - ~FullMask(:,:,1)*.1)
                colormap('jet')
                caxis([-.1 .4])
                
                set(gca,'LooseInset',get(gca,'TightInset'));
                saveMatchesName=['NMIPlot' rater '_f' num2str(f) '_m' num2str(m) '.bmp'];
                saveas(figID1,fullfile(outDirBase,saveMatchesName))
            end
        end
        
    end
end
%save results
save('OverlapSimilarityResults.mat','avgCorr_all',...
    'avgNMI_all','stdDevCorr_all',...
    'stdDevNMI_all','-v7.3');
