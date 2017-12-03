function  [outNameList, TotalTransform, MatchedTo_Followup_to_Baseline] = AOAlignToExistingMontage(imageDir_Baseline, posFileLoc_Baseline, baselineMontageSaved, imageDir_Followup, posFileLoc_Followup, outputDir, device_mode, ModalitiesSrchStrings,TransType,featureType)
%function for matching a followup montage onto a pre-built baseline montage
%Inputs
%imageDir_Baseline -- Image directory path for baseline images
%posFileLoc_Baseline -- Filepath for the nominal position file for the baseline images
%baselineMontageSaved -- Filepath for the 'AOMontageSave.mat' file from the completed baseline montage
%imageDir_Followup -- Image directory path for followup images
%posFileLoc_Followup -- Filepath for the nominal position file for the followup images
%outputDir -- Folder location for where the montaged images will be saved
%device_mode -- The device type used to acquire the images ['multi_modal' or 'canon']
%ModalitiesSrchStrings -- Cell array of search strings to find each modality
%       Example -- ModalitiesSrchStrings = {'confocal'; 'split_det'; 'avg'};
%TransType -- Index for the type of transformation used by the matching:
%  0 - Translation Only
%  1 - Translation + Rotation
%  2 - Rotation + Translation + Scaling
%  3 - Affine
%  4 - Rotation + Translation + Individual Scaling
%featureType -- Feature type to use for matching
%  0 - SIFT
%  1 - Constellation Feature

%Outputs:
%outNameList -- List of all file locations for saved montaged images
%TotalTransform -- Total transformation matrix for all followup images
%MatchedTo_Followup_to_Baseline -- list with index references from followup to baseline (which image in baseline each followup image matched to)

%Written by Min Chen (minchen1@upenn.edu)

tic
%Algorithm Parameters
NumOkMatchesThresh = 5; %Threshold for number of SIFT matches needed before accepting the transformation
matchexp = '_\d\d\d\d_ref_\d'; %String match that is expected to show up in the filename of each image. E.g. '_0018_ref_7_'
parallelFlag = exist('parfor');

%load baseline and followup data

[inData_Baseline, MN] = organizeDataByModality(imageDir_Baseline, ModalitiesSrchStrings);
[inData_Followup, MN] = organizeDataByModality(imageDir_Followup, ModalitiesSrchStrings);

MN = 1;

%number of images in each
N_Baseline = size(inData_Baseline,2);
N_Followup = size(inData_Followup,2);

%initialize all variables
RelativeTransformToRef_Followup = zeros(3,3,N_Followup);%stores relative transform between matched pairs
Matched_Followup = zeros(1,N_Followup);%store if image is matched already
MatchedTo_Followup_to_Baseline = zeros(1,N_Followup);%keeps track of which image in baseline each follow_up image is matched to

%stores all pairwise transforms
ResultsNumOkMatches = -ones(N_Followup,N_Baseline);
ResultsNumMatches = -ones(N_Followup,N_Baseline);
ResultsScaleToRef = -ones(N_Followup,N_Baseline);
ResultsTransformToRef = zeros(3,3,N_Followup,N_Baseline);

%read position file
[imageFilename_Baseline, eyeSide, pixelScale_Baseline, LocXY_Baseline, ID_Baseline, NC, errorFlag] = readPositionFile(imageDir_Baseline, inData_Baseline, posFileLoc_Baseline, device_mode, matchexp, MN, N_Baseline);
[imageFilename_Followup, eyeSide, pixelScale_Followup, LocXY_Followup, ID_Followup, NC, errorFlag] = readPositionFile(imageDir_Followup, inData_Followup, posFileLoc_Followup, device_mode, matchexp, MN, N_Followup);

%catch errors
if(errorFlag)
    errordlg(errorFlag);
    outNameList = [];
    return
end

%sort using LocXY
[LocXY_Baseline, inData_Baseline, imageFilename_Baseline, pixelScale_Baseline,ID_Baseline] = sortUsingLocXY(LocXY_Baseline, inData_Baseline, imageFilename_Baseline, pixelScale_Baseline, ID_Baseline, MN);
[LocXY_Followup, inData_Followup, imageFilename_Followup, pixelScale_Followup,ID_Followup] = sortUsingLocXY(LocXY_Followup, inData_Followup, imageFilename_Followup, pixelScale_Followup, ID_Baseline, MN);

%If appending we start by loading variable from the save
saved = load(baselineMontageSaved,'N','inData','TotalTransform','LocXY');
%Check saved info for baseline matches with loaded info for baseline
verified = 1;%verify if image was available in previous run
for n = 1:N_Baseline
    for m=1:MN
        if(~isequal(inData_Baseline(m,n),saved.inData(m,n)))% check if the name is same from previous results.
            verified = 0;
        end
    end
end
if(~verified)
    errordlg(['Baseline data do not previously built baseline montage']);
    outNameList = [];
    return
    
end


%calculate features
[f_all_Baseline, d_all_Baseline, h] = calculateFeatures(imageFilename_Baseline, parallelFlag, pixelScale_Baseline, featureType, MN, N_Baseline);
[f_all_Followup, d_all_Followup, h] = calculateFeatures(imageFilename_Followup, parallelFlag, pixelScale_Followup, featureType, MN, N_Followup);

%Begin Matching
waitbar(sum(Matched_Followup)/N_Followup,h,strcat('Aligning Images (',num2str(100*sum(Matched_Followup)/N_Followup,3),'%)'));

%find the closest unmatched image to the origin. Set as new reference frame
%    unMatchedI = find(Matched == 0);
%    unMatchedLocXY = LocXY(:,unMatchedI);
%    [sortUnMatchedLocXY, I] = sortrows(abs(unMatchedLocXY)',[1,2]);

%set reference frame
%    refIndex = unMatchedI(I(1));
%    RelativeTransformToRef(:,:,refIndex) = eye(3,3);
%    Matched(refIndex) = 1;
%    MatchedTo(refIndex) = refIndex;


%Allocate/reset Variables
bestRefIndex = 0;
bestTransform = zeros(3,3);
bestNumOkMatches = 0;
bestNumMatches = 1000000;
searchLevel = 1;
percentageThresh = .50;
stuckFlag = 0;

%find all matches linking to this reference frame
while(stuckFlag == 0)
    stuckFlag = 1; %make sure something changes this latest iteration
    for n=1:N_Followup
        if(Matched_Followup(n) == 0) % look for an unmatched Images
            disp(strcat('Matching n',int2str(n), ' at location (', num2str(LocXY_Followup(1,n)),',', ...
                num2str(LocXY_Followup(2,n)),')'));
            %reset variables
            bestNumOkMatches = 0;%Set to Min
            bestNumMatches = 1000000;%Set to Max
            
            for refIndex=1:N_Baseline % look at all possible references
                CurrentDistToRef = sum(abs(LocXY_Baseline(:,refIndex) - LocXY_Followup(:,n)));
                if(CurrentDistToRef <= searchLevel && refIndex ~= n)%Check reference is within one distance to target
                    
                    disp(strcat('  Checking against n',int2str(refIndex), ' at location (', num2str(LocXY_Baseline(1,refIndex)),',', ...
                        num2str(LocXY_Baseline(2,refIndex)),')'));
                    
                    if (ResultsNumOkMatches(n,refIndex) == -1) %check if match has already been checked before
                        refImg = cell(MN,1);
                        currentImg= cell(MN,1);
                        for m= 1:MN
                            refImg{m} = imresize(imread(char(imageFilename_Baseline{m,refIndex})), pixelScale_Baseline(refIndex),'bilinear' );
                            currentImg{m} = imresize( imread(char(imageFilename_Followup{m,n})), pixelScale_Followup(n),'bilinear' );
                        end
                        saveMatchesName = strcat('Matches_n',num2str(n),'_(', num2str(LocXY_Followup(1,n)),',', ...
                            num2str(LocXY_Followup(2,n)),')','_to_','n',num2str(refIndex),'_(', num2str(LocXY_Baseline(1,refIndex)),...
                            ',',num2str(LocXY_Baseline(2,refIndex)),')');
                        
                        saveMatchesName=fullfile(outputDir,saveMatchesName);

                        [relativeTransform, numOkMatches, numMatches, bestScale]=sift_mosaic_fast_MultiModal(refImg, currentImg, saveMatchesName,0,f_all_Baseline(:,refIndex),d_all_Baseline(:,refIndex),f_all_Followup(:,n),d_all_Followup(:,n),TransType,[],featureType);
                        
                        ResultsNumOkMatches(n,refIndex) = numOkMatches;
                        ResultsNumMatches(n,refIndex) = numMatches;
                        ResultsTransformToRef(:,:,n,refIndex) = relativeTransform;
                        ResultsScaleToRef(n,refIndex) =  bestScale;
                    end
                    
                    %set as best if better than current best
                    if(ResultsNumOkMatches(n,refIndex)>bestNumOkMatches || ((ResultsNumOkMatches(n,refIndex)==bestNumOkMatches) && ...
                            (ResultsNumOkMatches(n,refIndex)/ResultsNumMatches(n,refIndex) > bestNumOkMatches/bestNumMatches)))
                        bestRefIndex = refIndex;
                        bestTransform(:,:) = ResultsTransformToRef(:,:,n,refIndex);
                        bestNumOkMatches = ResultsNumOkMatches(n,refIndex);
                        bestNumMatches = ResultsNumMatches(n,refIndex);
                    end
                end
            end
            
            %use best match if it provides enough viable matches or there is no
            %more incomplete neighbors. Otherwise we go another round
            if(bestNumOkMatches >= NumOkMatchesThresh )%|| ((bestNumOkMatches/bestNumMatches) > percentageThresh && bestNumOkMatches > 1))
                disp(strcat('  Best match found n',int2str(bestRefIndex), ' at location (', num2str(LocXY_Baseline(1,bestRefIndex)),',', ...
                    num2str(LocXY_Baseline(2,bestRefIndex)),')'));
                RelativeTransformToRef_Followup(:,:,n) = bestTransform;
                Matched_Followup(n) = 1;
                MatchedTo_Followup_to_Baseline(n) = bestRefIndex;
                
                %save best one
                refIndex = bestRefIndex;
                saveMatchesName = strcat('Matches_n',num2str(n),'_(', num2str(LocXY_Followup(1,n)),',', ...
                            num2str(LocXY_Followup(2,n)),')','_to_','n',num2str(refIndex),'_(', num2str(LocXY_Baseline(1,refIndex)),...
                            ',',num2str(LocXY_Baseline(2,refIndex)),')');
                 saveMatchesName=fullfile(outputDir,saveMatchesName);
                 mkdir(saveMatchesName);

                 sift_mosaic_fast_MultiModal(refImg, currentImg, saveMatchesName,1,f_all_Baseline(:,refIndex),d_all_Baseline(:,refIndex),f_all_Followup(:,n),d_all_Followup(:,n),TransType,[],featureType);
                
                stuckFlag = 0;
                waitbar(sum(Matched_Followup)/N_Followup,h,strcat('Aligning Images (',num2str(100*sum(Matched_Followup)/N_Followup,3),'%)'));
            end
        end
    end
    
    if(stuckFlag == 1)%Release some constraints if stuck
        if(searchLevel == 1)%raise search radius
            searchLevel = 2;
            stuckFlag = 0;
%        elseif(percentageThresh > .25)%lower require match percentage
%            percentageThresh = percentageThresh - .1;
%            stuckFlag = 0;
        end
        %if both constraints have already been lowered, then we just start a new
        %branch
    end
end
waitbar(sum(Matched_Followup)/N_Followup,h,strcat('Alignment Complete (100%). Writing Outputs.'));

%find reference frames for different montage pieces

NumOfRefs = 1;

%Calc Total Transform relative to baseline Total Transform
TotalTransform = zeros(size(RelativeTransformToRef_Followup));
for n = 1:N_Followup
    if(Matched_Followup(n) == 1)
        H_Followup = RelativeTransformToRef_Followup(:,:,n);
        baselineRefIndex = MatchedTo_Followup_to_Baseline(n);
        H_Baseline =saved.TotalTransform(:,:,baselineRefIndex);
        TotalTransform(:,:,n) = H_Followup*H_Baseline;
    end
end

%calculate global bounding box
maxXAll =-1000000000;
minXAll =1000000000;
maxYAll =-1000000000;
minYAll =1000000000;

for n = 1:N_Followup
    if(Matched_Followup(n) == 1)
        
        im = imresize( imread(char(imageFilename_Followup{1,n})), pixelScale_Followup(n),'nearest'); % We don't care what they look like
        H = TotalTransform(:,:,n);
        
        %transform the 4 corners
        box = [1  size(im,2) size(im,2)  1 ;
            1  1           size(im,1)  size(im,1) ;
            1  1           1            1 ] ;
        box_ = inv(H) * box ;
        box_(1,:) = box_(1,:) ./ box_(3,:) ;
        box_(2,:) = box_(2,:) ./ box_(3,:) ;
        
        maxXAll = max([maxXAll box_(1,:)]);
        minXAll = min([minXAll box_(1,:)]);
        maxYAll = max([maxYAll box_(2,:)]);
        minYAll = min([minYAll box_(2,:)]);
    end
end

for n = 1:N_Baseline
    im = imresize( imread(char(imageFilename_Baseline{1,n})), pixelScale_Baseline(n),'nearest'); % We don't care what they look like
    H = saved.TotalTransform(:,:,n);
    
    %transform the 4 corners
    box = [1  size(im,2) size(im,2)  1 ;
        1  1           size(im,1)  size(im,1) ;
        1  1           1            1 ] ;
    box_ = inv(H) * box ;
    box_(1,:) = box_(1,:) ./ box_(3,:) ;
    box_(2,:) = box_(2,:) ./ box_(3,:) ;
    
    maxXAll = max([maxXAll box_(1,:)]);
    minXAll = min([minXAll box_(1,:)]);
    maxYAll = max([maxYAll box_(2,:)]);
    minYAll = min([minYAll box_(2,:)]);
end

%this is the image grid to output
ur = minXAll:maxXAll;
vr = minYAll:maxYAll;
[u,v] = meshgrid(ur,vr) ;

if(NumOfRefs == 1)
    outNameList = cell(1,MN);%Just 1 output for each modality if only one piece
else
    outNameList = cell(NumOfRefs+1,MN);%Otherwise one for each piece and extra one for all the pieces combined
end
numWritten=0;%keeps track of how many images written to disk


%save variables
%save(fullfile(outputDir,'AOMontageSave'),'LocXY','inData','TransType',...
%    'ResultsNumOkMatches','ResultsNumMatches',...
%    'ResultsTransformToRef','f_all','d_all','N_Followup','ResultsScaleToRef','MatchedTo','TotalTransform',...
%    'RelativeTransformToRef');

outputDirMatched = fullfile(outputDir,'matched');
outputDirUnmatched = fullfile(outputDir,'unmatched');

mkdir(outputDirMatched)
mkdir(outputDirUnmatched)

for m = 1:MN
    %initialize blank combined image of all pieces for the modality
    im =  imread(char(imageFilename_Followup{m,1}));
    imCombinedAll = vl_imwbackward(im2double(im),u,v);
    imCombinedAll=imCombinedAll(:,:,1);
    imCombinedAll(:,:,:) = 0;
    
    for i = 1: NumOfRefs
        %initialize blank combined image for the modality/piece
        im =  imread(char(imageFilename_Followup{m,1}));
        imCombined = vl_imwbackward(im2double(im),u,v);
        imCombined = imCombined(:,:,1);
        imCombined(:,:,:) = 0;
        
        for n = 1:N_Followup
            if(Matched_Followup(n) == 1)
                %read each image, and then transform
                im = imread(char(imageFilename_Followup{m,n}));
                H = TotalTransform(:,:,n);
                z_ = H(3,1) * u + H(3,2) * v + H(3,3);
                u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
                v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
                im_ = vl_imwbackward(im2double(im),u_,v_) ;
                im_ = im_(:,:,1);
                
                %add to combined image
                nonzero = im_>0;
                imCombined(nonzero) = im_(nonzero);
                %imCombined=max(imCombined, im_);
                %figure(i)
                %imshow(imCombined)
                
                %save each individually transformed image
                [pathstr,name,ext] = fileparts(char(imageFilename_Followup{m,n})) ;
                
                rgba = repmat(im_(:,:,1),[1,1,2]);
                rgba(:,:,2) = ~isnan(im_(:,:,1));
                rgba = uint8(round(rgba*255));
                
                saveFileName=[name,'_aligned_to_ref',num2str(i),'_m',num2str(m),'.tif'];
                saveTif(rgba,outputDirMatched,saveFileName);
            else
                copyfile(imageFilename_Followup{m,n},outputDirUnmatched);
            end
            numWritten = numWritten+1;
            waitbar(numWritten/(N_Followup*MN),h,strcat('Writing Outputs (',num2str(100*numWritten/(N_Followup*MN),3),'%)'));
        end
        
        %add to all combined image
        nonzero = imCombined>0;
        imCombinedAll(nonzero) = imCombined(nonzero);
        
        
        %figure(i);clf;
        %imshow(imCombined)
        %save combined image for each piece
        if(NumOfRefs > 1)%only necessary if more than one piece
            saveFileName = strcat('ref_',num2str(i),'_combined_m',num2str(m),'.tif');
            
           % if(AppendToExisting)
           %     outNameList{i+1,m}=fullfile('Append',saveFileName);
           % else
                outNameList{i+1,m}=saveFileName;
           % end
            
            rgba = repmat(imCombined(:,:,1),[1,1,2]);
            rgba(:,:,2) = ~isnan(imCombined(:,:,1));
            rgba = uint8(round(rgba*255));
            saveTif(rgba,outputDir,saveFileName);
        end
    end
    
    %save the combined image of all the pieces
    saveFileName = strcat('all_ref_combined_m',num2str(m),'.tif');
  %  if(AppendToExisting)
  %      outNameList{1,m}=fullfile('Append',saveFileName);
  %  else
        outNameList{1,m}=saveFileName;
  %  end
    
    rgba = repmat(imCombinedAll(:,:,1),[1,1,2]);
    rgba(:,:,2) = ~isnan(imCombinedAll(:,:,1));
    rgba = uint8(round(rgba*255));
    saveTif(rgba,outputDir,saveFileName);
    
end
runtime=toc


