function  outNameList = AOMosiacAllMultiModal(imageDir, posFileLoc, outputDir, device_mode, ModalitiesSrchStrings,TransType,AppendToExisting,MontageSave)
%Main AO Montaging Function that creates a full montage from an input
%directory with images and nominal coordinate location
%Inputs:
%imageDir -- Folder location containing all input images to be montaged 
%posFileLoc -- An excel spreadsheet file indicating the coordinate
%              locations of each image. [See demo->AOIP for example file.]
%outputDir -- Folder location for where the montaged images will be saved
%device_mode -- The device type used to acquire the images [aoip or canon] 
%ModalitiesSrchStrings -- Cell array of search strings to find each modality 
%       Example -- ModalitiesSrchStrings = {'confocal'; 'split_det'; 'avg'};
%TransType -- Index for the type of transformation used by the matching:
%  0 - Translation Only
%  1 - Translation + Rotation
%  2 - Rotation + Translation + Scaling
%  3 - Affine
%  4 - Rotation + Translation + Individual Scaling
%AppendToExisting -- Flag indicating if the montage is appending an
%                    existing montage [0 or 1]
%MontageSave -- File location of existing AOMontageSave.m to be appended to.
%               Only applicable if AppendToExisting = 1;

%Outputs:
%outNameList -- List of all file locations for saved montaged images
%
%
%Written by Min Chen (minchen1@upenn.edu)
%Edits by Robert F. Cooper (rfcooper@sas.upenn.edu)

tic
%Load Data
%baseDir = 'C:\Users\Min\Dropbox\AOSLOImageProcessing\ConstructMontage\InputImages_Set1\';
%fileID = fopen(fullfile(baseDir,'ImageInfo.txt'));

%Load info from descriptor file
% formatSpec = '%s';
% C_text = textscan(fileID,formatSpec,9,'Delimiter','\t');
% C = textscan(fileID, '%d8 %s %s %s %s %d8 %s %d8 %d8', 'delimiter', '\t');
%N = size(C{1},1);
parallelFlag = exist('parfor');
Allfiles = dir(fullfile(imageDir,'*.tif'));
Allfiles = {Allfiles.name};

MN = size(ModalitiesSrchStrings,1);
inData = [];

%create multi-modal data structure
%depending on parameter array modalitiesToUse
counter = 0;
for m = 1:MN
    
    if (~isempty(ModalitiesSrchStrings{m}))%check it's not empty
        imagesFound = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, ModalitiesSrchStrings{m}))));%search for image of this modality
        if(~isempty(imagesFound))
            inData = [inData; imagesFound];
            counter = counter + 1;
        end

    end
    %confocal = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_confocal_'))));
    %darkfield = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_avg_'))));
    %splitDecision = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_split_'))));
    
end
MN = counter;%only using nonempty identifiers
N = size(inData,2);


%initialize all variables
LocXY = nan(2,N);%set NaN to start
RelativeTransformToRef = zeros(3,3,N);%stores relative transform between matched pairs
imageFilename = cell(MN,N);%stores all filenames
Matched = zeros(1,N);%store if image is matched already
MatchedTo = zeros(1,N);%keeps track of which image is matched to which

%stores all pairwise transforms
ResultsNumOkMatches = -ones(N,N);
ResultsNumMatches = -ones(N,N);
ResultsTransformToRef = zeros(3,3,N,N);

%stores sift features for each image
f_all = cell(MN,N);
d_all = cell(MN,N);

if strcmp(device_mode, 'aoip')

    %load position info from excel spreadsheet
    [temp,temp,C] = xlsread(posFileLoc);


    %verify that the image id's line up for all modalities
    %example _0018_ref_7_
    matchexp = '_\d\d\d\d_ref_\d';
    eyeSide = 'OS';
    for n = 1:N
        %build filename structure and check to make sure all modalities are
        %present for all ids
        imageFilename{1,n} = fullfile(imageDir, inData{1,n});
        ImageID_m1 = regexpi(inData{1,n},matchexp,'match');
        for m = 2:MN
            imageFilename{m,n} = fullfile(imageDir, inData{m,n});
            ImageID_mf = regexpi(inData{m,n},matchexp,'match');
            if(~strcmpi(ImageID_m1,ImageID_mf));
                errordlg(['Error: Mismatch detected. Every image number must have the same number of modalities. Check image ' ImageID_m1]);
                outNameList = [];
                return
            end
        end
        %match with info from excel
        for i = 1:size(C,1)

            if(strcmpi(C{i,1},'eye'))
                eyeSide = C{i,2};
            end

            if (~isempty(strfind(inData{1,n}, C{i,1})))

                %first try looking at coordinate grid
                if(size(C,2) >= 3)
                    Loc = strsplit(C{i,3},',');
                    if(size(Loc,2) == 2)
                        LocXY(1,n) = str2double(strtrim(Loc{1}));
                        LocXY(2,n) = str2double(strtrim(Loc{2}));
                    end
                end

                if(size(C,2) >= 2)
                    %coordinate grind c
                    if(isnan(LocXY(1,n)) || isnan(LocXY(2,n))) 
                        LocXY(:,n) = parseShorthandLoc(C{i,2},eyeSide);
                    end
                end
                if(isnan(LocXY(1,n)) || isnan(LocXY(2,n)))
                    errordlg(['Error: Location missing or invalid for image ' ImageID_m1]);
                    outNameList = [];
                    return
                end
                break;
            end
        end

    end
    
elseif strcmp(device_mode, 'canon')
    for n = 1:N
        [ eyeSide, LocXY(:,n) ] = parseCanonFName( inData{1,n} );
        imageFilename{1,n} = fullfile(imageDir, inData{1,n});
    end
end

%sort using LocXY
[sortLocXY I] = sortrows(abs(LocXY)',[1,2]);
LocXY=LocXY(:,I);

for m = 1:MN
    imageFilename(m,:)=imageFilename(m,I);
end

if(~AppendToExisting)
    %If new then calculate All Sift Features
    h = waitbar(0,'Calculating Sift Features (0%)');
    for n=1:N
        if(parallelFlag)
            parfor m = 1:MN
                im = im2single(imread(char(imageFilename{m,n})));
                [f1,d1] = vl_sift(im,'Levels',55);
                f_all{m,n} = f1;
                d_all{m,n} = d1;
            end
        else
            for m = 1:MN
                im = im2single(imread(char(imageFilename{m,n})));
                [f1,d1] = vl_sift(im,'Levels',55);
                f_all{m,n} = f1;
                d_all{m,n} = d1;
            end
        end
        waitbar(n/(N),h,strcat('Calculating Sift Features (',num2str(100*n/N,3),'%)'));
    end
else
    %If appending we start by loading variable from the save
    saved = load(MontageSave,'LocXY','inData','TransType',...
        'ResultsNumOkMatches','ResultsNumMatches',...
        'ResultsTransformToRef','f_all','d_all','N');
    if(~isfield(saved,'TransType') || (saved.TransType ~= TransType))
        choice = questdlg('Selected transformation type does not match that saved in AOMontageSave.mat. Proceed anyways?','Warrning');
        if(~isequal(choice,'Yes'))
            outNameList = [];
            return;
        end
    end
    
    outputDir = fullfile(outputDir,'Append');
    mkdir(outputDir);
    
    %Check saved info and match indexes
    verified = zeros(1,N);%verify if image was available in previous run
    transferFrom = zeros(1,N);%store reference index from previous run of same image
    for n = 1:N
        for s = 1:saved.N
            for m=1:MN
                if(isequal(inData(1,n),saved.inData(m,s)) && isequal(LocXY(:,n),saved.LocXY(:,s)))
                    verified(n) = 1;
                    transferFrom(n) = s;
                end
            end
            
        end
    end
    
    L = length(find(verified~=1));
    newCounter=0;
    %transfer saved info and calculate sift features for new images
    h = waitbar(0,'Calculating Sift Features For New Images (0%)');
    for n1 = 1:N
        if(verified(n1))%check if it exists in save
            s1 = transferFrom(n1);
            for m = 1:MN%copy sift features
                f_all{m,n1} = saved.f_all{m,s1};
                d_all{m,n1} = saved.d_all{m,s1};
            end
            for n2 = 1:N%copy pairwise matches
                if(verified(n2))
                    s2 = transferFrom(n2);
                    ResultsNumOkMatches(n1,n2) = saved.ResultsNumOkMatches(s1,s2);
                    ResultsNumMatches(n1,n2) = saved.ResultsNumMatches(s1,s2);
                    ResultsTransformToRef(:,:,n1,n2) = saved.ResultsTransformToRef(:,:,s1,s2);
                end
            end
        else%if not calculate sift features
            if(parallelFlag)
                parfor m = 1:MN
                    im = im2single(imread(char(imageFilename{m,n1})));
                    [f1,d1] = vl_sift(im,'Levels',55);
                    f_all{m,n1} = f1;
                    d_all{m,n1} = d1;
                end
            else
                for m = 1:MN
                    im = im2single(imread(char(imageFilename{m,n1})));
                    [f1,d1] = vl_sift(im,'Levels',55);
                    f_all{m,n1} = f1;
                    d_all{m,n1} = d1;
                end
            end
            
            waitbar(newCounter/(L),h,strcat('Calculating Sift Features For New Images(',num2str(100*newCounter/L,3),'%)'));
        end
    end
    
end

%Begin Matching
while (sum(Matched) < N)
    
    waitbar(sum(Matched)/N,h,strcat('Aligning Images (',num2str(100*sum(Matched)/N,3),'%)'));
    
    %find the closest unmatched image to the origin. Set as new reference frame
    unMatchedI = find(Matched == 0);
    unMatchedLocXY = LocXY(:,unMatchedI);
    [sortUnMatchedLocXY, I] = sortrows(abs(unMatchedLocXY)',[1,2]);
    
    %set reference frame
    refIndex = unMatchedI(I(1));
    RelativeTransformToRef(:,:,refIndex) = eye(3,3);
    Matched(refIndex) = 1;
    MatchedTo(refIndex) = refIndex;
    
    
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
        for n=1:N
            if(Matched(n) == 0) % look for an unmatched Images
                disp(strcat('Matching n',int2str(n), ' at location (', num2str(LocXY(1,n)),',', ...
                    num2str(LocXY(2,n)),')'));
                %reset variables
                bestNumOkMatches = 0;%Set to Min
                bestNumMatches = 1000000;%Set to Max
                for refIndex=1:N % look at all possible references
                    CurrentDistToRef = sum(abs(LocXY(:,refIndex) - LocXY(:,n)));
                    if(CurrentDistToRef <= searchLevel && refIndex ~= n)%Check reference is within one distance to target
                        if(Matched(refIndex) == 1)%Check that reference has been matched
                            
                            disp(strcat('  Checking against n',int2str(refIndex), ' at location (', num2str(LocXY(1,refIndex)),',', ...
                                num2str(LocXY(2,refIndex)),')'));
                            
                            if (ResultsNumOkMatches(n,refIndex) == -1) %check if match has already been checked before
                                refImg = cell(MN,1);
                                currentImg= cell(MN,1);
                                for m= 1:MN
                                    refImg{m} = imread(char(imageFilename{m,refIndex}));
                                    currentImg{m} = imread(char(imageFilename{m,n}));
                                end
                                saveMatchesName = strcat('Matches_n',num2str(n),'_(', num2str(LocXY(1,n)),',', ...
                                    num2str(LocXY(2,n)),')','_to_','n',num2str(refIndex),'_(', num2str(LocXY(1,refIndex)),...
                                    ',',num2str(LocXY(2,refIndex)),')','.jpg');
                                
                                saveMatchesName=fullfile(outputDir,saveMatchesName);
                                [relativeTransform, numOkMatches, numMatches]=sift_mosaic_fast_MultiModal(refImg, currentImg, saveMatchesName,0,f_all(:,refIndex),d_all(:,refIndex),f_all(:,n),d_all(:,n),TransType);
                                
                                ResultsNumOkMatches(n,refIndex) = numOkMatches;
                                ResultsNumMatches(n,refIndex) = numMatches;
                                ResultsTransformToRef(:,:,n,refIndex) = relativeTransform;
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
                end
                
                %use best match if it provides enough viable matches or there is no
                %more incomplete neighbors. Otherwise we go another round
                if(bestNumOkMatches > 10 || ((bestNumOkMatches/bestNumMatches) > percentageThresh && bestNumOkMatches > 1))
                    disp(strcat('  Best match found n',int2str(bestRefIndex), ' at location (', num2str(LocXY(1,bestRefIndex)),',', ...
                        num2str(LocXY(2,bestRefIndex)),')'));
                    RelativeTransformToRef(:,:,n) = bestTransform;
                    Matched(n) = 1;
                    MatchedTo(n) = bestRefIndex;
                    stuckFlag = 0;
                    waitbar(sum(Matched)/N,h,strcat('Aligning Images (',num2str(100*sum(Matched)/N,3),'%)'));
                end
            end
        end
        
        if(stuckFlag == 1)%Release some constraints if stuck
            if(searchLevel == 1)%raise search radius
                searchLevel = 2;
                stuckFlag = 0;
            elseif(percentageThresh > .25)%lower require match percentage
                percentageThresh = percentageThresh - .1;
                stuckFlag = 0;
            end
            %if both constraints have already been lowered, then we just start a new
            %branch
        end
    end
end


AllRefIndex = [];

for n = 1:N
    if(MatchedTo(n) == n)
        AllRefIndex = [AllRefIndex n];
    end
end

NumOfRefs = size(AllRefIndex,2);
RefChains = cell(NumOfRefs,1);

%Calc Total Transform
TotalTransform = zeros(size(RelativeTransformToRef));
for n = 1:N
    H = RelativeTransformToRef(:,:,n);
    nextRef = MatchedTo(n);
    
    %Find Total Transform
    while (nextRef ~= MatchedTo(nextRef))
        nextH = RelativeTransformToRef(:,:,nextRef);
        H = H*nextH;
        nextRef = MatchedTo(nextRef);
    end
    
    %Add to refence chain
    ind = find(nextRef == AllRefIndex,1);
    RefChains{ind} = [RefChains{ind} n];
    
    TotalTransform(:,:,n) = H;
end

outNameList = cell(NumOfRefs,MN);

for m = 1:MN
    for i = 1: NumOfRefs
        %find max and min bounding box over all iamges attached to this ref
        
        maxX =-1000000000;
        minX =1000000000;
        maxY =-1000000000;
        minY =1000000000;
        
        for n = 1:N
            im = imread(char(imageFilename{m,n}));
            H = TotalTransform(:,:,n);
            
            %transform the 4 corners
            box = [1  size(im,2) size(im,2)  1 ;
                1  1           size(im,1)  size(im,1) ;
                1  1           1            1 ] ;
            box_ = inv(H) * box ;
            box_(1,:) = box_(1,:) ./ box_(3,:) ;
            box_(2,:) = box_(2,:) ./ box_(3,:) ;
            
            maxX = max([maxX box_(1,:)]);
            minX = min([minX box_(1,:)]);
            maxY = max([maxY box_(2,:)]);
            minY = min([minY box_(2,:)]);
        end
        
        ur = minX:maxX;
        vr = minY:maxY;
        [u,v] = meshgrid(ur,vr) ;
        
        im =  imread(char(imageFilename{m,AllRefIndex(i)}));
        imCombined = vl_imwbackward(im2double(im),u,v);
        
        for n = RefChains{i}
            im = imread(char(imageFilename{m,n}));
            H = TotalTransform(:,:,n);
            z_ = H(3,1) * u + H(3,2) * v + H(3,3);
            u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
            v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
            im_ = vl_imwbackward(im2double(im),u_,v_) ;
            
            nonzero = im_>0;
            imCombined(nonzero) = im_(nonzero);
            %imCombined=max(imCombined, im_);
            %figure(i)
            %imshow(imCombined)
            [pathstr,name,ext] = fileparts(char(imageFilename{m,n})) ;
            
            rgba = repmat(im_,[1,1,2]);
            rgba(:,:,2) = ~isnan(im_);
            rgba = uint8(round(rgba*255));
            
            %# create a tiff object
            tob = Tiff(fullfile(outputDir,[name,'_aligned_to_ref',num2str(i),'_m',num2str(m),'.tif']),'w');
            
            %# you need to set Photometric before Compression
            tob.setTag('Photometric',Tiff.Photometric.MinIsBlack)
            tob.setTag('Compression',Tiff.Compression.LZW)
            
            %# tell the program that channel 4 is alpha
            tob.setTag('ExtraSamples',Tiff.ExtraSamples.AssociatedAlpha)
            
            %# set additional tags (you may want to use the structure
            %# version of this for convenience)
            tob.setTag('ImageLength',size(im_,1));
            tob.setTag('ImageWidth',size(im_,2));
            tob.setTag('BitsPerSample',8);
            tob.setTag('RowsPerStrip',16);
            tob.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Separate);
            tob.setTag('Software','MATLAB')
            tob.setTag('SamplesPerPixel',2);
            
            %# write and close the file
            tob.write(rgba)
            tob.close
            
        end
        
        %figure(i);clf;
        %imshow(imCombined)
        outName = strcat('ref_',num2str(i),'_combined_m',num2str(m),'.tif');

        if(AppendToExisting)
            outNameList{i,m}=fullfile('Append',outName);
        else
            outNameList{i,m}=outName;
        end
        
        rgba = repmat(imCombined,[1,1,2]);
        rgba(:,:,2) = ~isnan(imCombined);
        rgba = uint8(round(rgba*255));
        
        %# create a tiff object
        tob = Tiff(fullfile(outputDir,outName),'w');
        
        %# you need to set Photometric before Compression
        tob.setTag('Photometric',Tiff.Photometric.MinIsBlack)
        tob.setTag('Compression',Tiff.Compression.LZW)
        
        %# tell the program that channel 4 is alpha
        tob.setTag('ExtraSamples',Tiff.ExtraSamples.AssociatedAlpha)
        
        %# set additional tags (you may want to use the structure
        %# version of this for convenience)
        tob.setTag('ImageLength',size(im_,1));
        tob.setTag('ImageWidth',size(im_,2));
        tob.setTag('BitsPerSample',8);
        tob.setTag('RowsPerStrip',16);
        tob.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Separate);
        tob.setTag('Software','MATLAB')
        tob.setTag('SamplesPerPixel',2);
        
        %# write and close the file
        tob.write(rgba)
        tob.close
        
    end
end

runtime=toc
%save variables
save(fullfile(outputDir,'AOMontageSave'),'LocXY','inData','TransType',...
    'ResultsNumOkMatches','ResultsNumMatches',...
    'ResultsTransformToRef','f_all','d_all','N');

