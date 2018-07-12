function  [outNameList, imageFilename, TotalTransform, f_all, d_all] = AOMosiacAllMultiModal(imageDir, posFileLoc, outputDir, device_mode, ModalitiesSrchStrings,TransType,AppendToExisting,MontageSave,export_to_pshop,featureType)
%Main AO Montaging Function that creates a full montage from an input
%directory with images and nominal coordinate location
%Inputs:
%imageDir -- Folder location containing all input images to be montaged
%posFileLoc -- An excel spreadsheet file indicating the coordinate
%              locations of each image. [See demo->Multi-Modal for example file.]
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
%AppendToExisting -- Flag indicating if the montage is appending an
%                    existing montage [0 or 1]
%MontageSave -- File location of existing AOMontageSave.m to be appended to.
%               Only applicable if AppendToExisting = 1;
%export_to_pshop -- flag for whether to also save as a photoshop montage
%featureType -- Feature type to use for matching
%  0 - SIFT
%  1 - Constellation Feature

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


%Algorithm Parameters
NumOkMatchesThresh = 10; %Threshold for number of SIFT matches needed before accepting the transformation
matchexp = '_\d\d\d\d_ref_\d'; %String match that is expected to show up in the filename of each image. E.g. '_0018_ref_7_'
if ~exist('featureType','var') || isempty(featureType)%default featuretype
    featureType=0;
end



%Load info from descriptor file
parallelFlag = exist('parfor');

%load data
[inData, MN] = organizeDataByModality(imageDir, ModalitiesSrchStrings);
N = size(inData,2);

%initialize all variables
RelativeTransformToRef = zeros(3,3,N);%stores relative transform between matched pairs
Matched = zeros(1,N);%store if image is matched already
MatchedTo = zeros(1,N);%keeps track of which image is matched to which

%stores all pairwise transforms
ResultsNumOkMatches = -ones(N,N);
ResultsNumMatches = -ones(N,N);
ResultsScaleToRef = -ones(N,N);
ResultsTransformToRef = zeros(3,3,N,N);

%read position file
[imageFilename, eyeSide, pixelScale, LocXY, ID, NC, errorFlag] = readPositionFile(imageDir, inData, posFileLoc, device_mode, matchexp, MN, N);
%catch errors
if(errorFlag)
    errordlg(errorFlag);
    outNameList = [];
    return
end

%sort using LocXY
[LocXY, inData, imageFilename, pixelScale,ID] = sortUsingLocXY(LocXY, inData, imageFilename, pixelScale, ID, MN);

if(~AppendToExisting)
    %If new then calculate All Features
    [f_all, d_all, h] = calculateFeatures(imageFilename, parallelFlag, pixelScale, featureType, MN, N);
else
    %If appending we start by loading variable from the save
    saved = load(MontageSave,'LocXY','inData','TransType',...
        'ResultsNumOkMatches','ResultsNumMatches',...
        'ResultsTransformToRef','f_all','d_all','N', 'ResultsScaleToRef');
        if(~isfield(saved,'TransType') || (saved.TransType ~= TransType))
            choice = questdlg('Selected transformation type does not match that saved in AOMontageSave.mat. Proceed anyways?','Warning');
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
                if(isequal(inData(1,n),saved.inData(m,s)))% check if the name is same from previous results.
                    LocXY(:,n) = saved.LocXY(:,s);
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
                    ResultsScaleToRef(n1,n2) = saved.ResultsScaleToRef(s1,s2);
                end
            end
        else%if not calculate sift features
            if(parallelFlag)
                
                parfor m = 1:MN
                    im = imread(char(imageFilename{m,n1}));
                    im = imresize( im2single( im(:,:,1) ), pixelScale(n),'bilinear' );
                    if(featureType == 0)
                        [f1,d1] = vl_sift(im,'Levels',SiftLevel);
                    elseif(featureType == 1)
                        [f1,d1] = gridFeatures(im(:,:,1));
                    else
                        [f1,d1] = vl_sift(im,'Levels',SiftLevel);
                    end
                    [f1_crop, d1_crop] = filterSiftFeaturesByROI(im, f1, d1, ROICropPct);
                    f_all{m,n1} = f1_crop;
                    d_all{m,n1} = d1_crop;
                end
                
            else
                for m = 1:MN
                    im = imread(char(imageFilename{m,n1}));
                    im = imresize( im2single( im(:,:,1) ), pixelScale(n),'bilinear' );
                    if(featureType == 0)
                        [f1,d1] = vl_sift(im,'Levels',SiftLevel);
                    elseif(featureType == 1)
                        [f1,d1] = gridFeatures(im(:,:,1));
                    else
                        [f1,d1] = vl_sift(im,'Levels',SiftLevel);
                    end
                    
                    [f1_crop, d1_crop] = filterSiftFeaturesByROI(im, f1, d1, ROICropPct);
                    f_all{m,n1} = f1_crop;
                    d_all{m,n1} = d1_crop;
                end
            end
            
            waitbar(newCounter/(L),h,strcat('Calculating Sift Features For New Images(',num2str(100*newCounter/L,3),'%)'));
        end
    end
    
end
if(parallelFlag)
    delete(gcp('nocreate'))
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
            if(Matched(n) == 0) && all(~isnan(LocXY(1,n))) % look for an unmatched Images, only match if they exist though...
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
                                    im = imread(char(imageFilename{m,refIndex}));
                                    refImg{m} = imresize( im(:,:,1), pixelScale(refIndex),'bilinear' );
                                    im = imread(char(imageFilename{m,n}));
                                    currentImg{m} = imresize( im(:,:,1), pixelScale(n),'bilinear' );
                                end
                                saveMatchesName = strcat('Matches_n',num2str(n),'_(', num2str(LocXY(1,n)),',', ...
                                    num2str(LocXY(2,n)),')','_to_','n',num2str(refIndex),'_(', num2str(LocXY(1,refIndex)),...
                                    ',',num2str(LocXY(2,refIndex)),')','.jpg');
                                
                                saveMatchesName=fullfile(outputDir,saveMatchesName);
                                [relativeTransform, numOkMatches, numMatches, bestScale]=sift_mosaic_fast_MultiModal(refImg, currentImg, saveMatchesName,0,f_all(:,refIndex),d_all(:,refIndex),f_all(:,n),d_all(:,n),TransType,[],featureType);
                                
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
                end
                
                %use best match if it provides enough viable matches or there is no
                %more incomplete neighbors. Otherwise we go another round
                if(bestNumOkMatches > NumOkMatchesThresh || ((bestNumOkMatches/bestNumMatches) > percentageThresh && bestNumOkMatches > 1))
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
waitbar(sum(Matched)/N,h,strcat('Alignment Complete (100%). Writing Outputs.'));

%find reference frames for different montage pieces
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

%translate discontinuous pieces of the montage relative to each other

%Centers of mass for each piece
CoMX = zeros(NumOfRefs,1);
CoMY = zeros(NumOfRefs,1);

%boundingbox for each piece
maxXRef = -1000000000*ones(NumOfRefs,1);
minXRef = 1000000000*ones(NumOfRefs,1);
maxYRef = -1000000000*ones(NumOfRefs,1);
minYRef = 1000000000*ones(NumOfRefs,1);
Xwidth = zeros(NumOfRefs,1);
Ywidth = zeros(NumOfRefs,1);

for i = 1: NumOfRefs
    
    %calculate centers of mass for each piece using nominal locations
    %We use median here to get cleaner locations
    CoMX(i) = median(LocXY(1,RefChains{i}));
    CoMY(i) = median(LocXY(2,RefChains{i}));
    
    %calculate bounding box for each piece
    for n = RefChains{i}
        if ~isempty(imageFilename{1,n})
            im = imread(char(imageFilename{1,n}));
            im = imresize( im(:,:,1), pixelScale(n),'nearest' ); % we don't care what they look like, only how big they are
            H = TotalTransform(:,:,n);

            %transform the 4 corners
            box = [1  size(im,2) size(im,2)  1 ;
                1  1           size(im,1)  size(im,1) ;
                1  1           1            1 ] ;
            box_ = inv(H) * box ;
            box_(1,:) = box_(1,:) ./ box_(3,:) ;
            box_(2,:) = box_(2,:) ./ box_(3,:) ;

            maxXRef(i) = max([maxXRef(i) box_(1,:)]);
            minXRef(i) = min([minXRef(i) box_(1,:)]);
            maxYRef(i) = max([maxYRef(i) box_(2,:)]);
            minYRef(i) = min([minYRef(i) box_(2,:)]);
        end
    end
    Xwidth(i) = maxXRef(i)-minXRef(i);
    Ywidth(i)= maxYRef(i)-minYRef(i);
end
%we want to make sure there is not two identical CoM, if there is, we push
%it in the larger direction.
ee = .0001;
checked=zeros(NumOfRefs,1);
for s = 1:NumOfRefs
    if(~checked(s))
        counter = 0;%keep tracks of how many duplicates, first one can be itself, so start at zero
        checked(s) = 1;
        for ss = 1:NumOfRefs
            if((CoMX(s) == CoMX(ss)) && (CoMY(s) == CoMY(ss)))
                checked(ss) = 1;
                if( CoMX(ss) >= CoMY(ss))%bump duplicate in the larger direction
                    CoMX(ss)=CoMX(ss)+ee*counter;
                else
                    CoMY(ss)=CoMY(ss)+ee*counter;
                end
                counter = counter + 1;%increment duplicate count
            end
        end
    end
    
end


%find relative translation between pieces using centers of mass
[refOrderX, refOrderX_I] = sort(CoMX);%find ascending order Index for each piece in X direction
[refOrderY, refOrderY_I] = sort(CoMY,'descend');%find descending order Index for each piece in Y direction


%find total translation for each piece
refGlobalTransX = zeros(NumOfRefs,1);
refGlobalTransY = zeros(NumOfRefs,1);

pad=30;
maxWidthX=Xwidth(refOrderX_I(1));
maxWidthY=Ywidth(refOrderY_I(1));

%if pieces have same nominal location in one direction, then we set the bounding box in the direction to min/max
%of all pieces at that location, this makes a more organized picture
%
for i = 1:NumOfRefs %ToDo: not very efficient... but small number of pieces in general
    ind = find(CoMX == CoMX(i));%find where location matches exactly
    maxXRef(ind) = max(maxXRef(ind));%set to max
    minXRef(ind) = min(minXRef(ind));%set to min
    
    ind = find(CoMY == CoMY(i));%find where location matches exactly
    maxYRef(ind) = max(maxYRef(ind));%set to max
    minYRef(ind) = min(minYRef(ind));%set to min
end



for s = 2:NumOfRefs %no need to translate the first one, start at 2
    
    
    %add width and pad for all previous pieces that do not share a number
    %The center of each piece starts at the same location (0,0), so the width we add is
    %max loc of the bounding of the previous image, and subtract min loc of the current image
    if(CoMX(refOrderX_I(s)) == CoMX(refOrderX_I(s-1)))%if same location as previous, use previous translation
        refGlobalTransX(refOrderX_I(s)) = refGlobalTransX(refOrderX_I(s-1));
        
    else%if different location, then use previous translation + max loc of previous location + pad - min loc of own
        refGlobalTransX(refOrderX_I(s)) = refGlobalTransX(refOrderX_I(s-1)) ...
            + maxXRef(refOrderX_I(s-1)) - minXRef(refOrderX_I(s)) + pad;
    end
    
    if(CoMY(refOrderY_I(s)) == CoMY(refOrderY_I(s-1)))
        refGlobalTransY(refOrderY_I(s)) = refGlobalTransY(refOrderY_I(s-1));
    else
        refGlobalTransY(refOrderY_I(s)) = refGlobalTransY(refOrderY_I(s-1)) ...
            + maxYRef(refOrderY_I(s-1)) - minYRef(refOrderY_I(s)) + pad;
    end
end

%Now adjust each transformation according to which piece they're in
for i = 1: NumOfRefs
    refGlobalTrans = eye(3,3);
    refGlobalTrans(1,3) = -refGlobalTransX(i); %these are pullback transforms, so negative of the distance you want to move
    refGlobalTrans(2,3) = -refGlobalTransY(i); %these are pullback transforms, so negative of the distance you want to move
    for n = RefChains{i}
        TotalTransform(:,:,n) = TotalTransform(:,:,n)*refGlobalTrans;
    end
end


%calculate global bounding box
maxXAll =-1000000000;
minXAll =1000000000;
maxYAll =-1000000000;
minYAll =1000000000;

for n = 1:N
    if ~isempty(imageFilename{1,n})
        im = imread(char(imageFilename{1,n}));
        im = imresize( im(:,:,1), pixelScale(n),'nearest' ); % we don't care what they look like, only how big they are        
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

%this is the image grid to output
ur = minXAll:maxXAll;
vr = minYAll:maxYAll;
% [u,v] = meshgrid(ur,vr) ;

if(NumOfRefs == 1)
    outNameList = cell(1,MN);%Just 1 output for each modality if only one piece
else
    outNameList = cell(NumOfRefs+1,MN);%Otherwise one for each piece and extra one for all the pieces combined
end
numWritten=0;%keeps track of how many images written to disk

fovlist = unique(pixelScale);
fovlist = cellfun(@(n) num2str(n,'%0.2f'), num2cell(round(fovlist,2)),'UniformOutput',false);

%save variables
save(fullfile(outputDir,'AOMontageSave'),'LocXY','inData','TransType',...
    'ResultsNumOkMatches','ResultsNumMatches',...
    'ResultsTransformToRef','f_all','d_all','N','ResultsScaleToRef','MatchedTo','TotalTransform',...
    'RelativeTransformToRef');

%% Determine the dominant direction of each shifted image
if export_to_pshop
    group_directions = cell(NumOfRefs,1);
    for i = 1: NumOfRefs
        for n = RefChains{i}
            
            % If we had 3 columns and we're not using the canon, then we can determine which should go in the
            % "fovea" bin.
            if ~strcmp(device_mode, 'canon') && (NC >= 3)
                if any(strcmpi(strtrim(ID{n}), {'TRC','TR','MTE','MT','TM','TLC','TL',...
                        'MRE','MR','C','CENTER','MLE','ML'...
                        'BRC','BR','MBE','MB','BM','BLC','BL',}))
                    group_directions{n} = 'Fovea';
                    continue;
                end
            end
            
            H = TotalTransform(:,:,n);
            
            [greaterdist ind] = max( abs(H(1:2,3)) );
            
            if ind==1
                if sign(H(ind,3)) == 1
                    if strcmpi(eyeSide,'os')
                        group_directions{n} = 'Temporal';
                    elseif strcmpi(eyeSide,'od')
                        group_directions{n} = 'Nasal';
                        
                    end
                else
                    if strcmpi(eyeSide,'os')
                        group_directions{n} = 'Nasal';
                    elseif strcmpi(eyeSide,'od')
                        group_directions{n} = 'Temporal';
                    end
                end
            else
                if sign(H(ind,3)) == 1
                    group_directions{n} = 'Superior';
                else
                    group_directions{n} = 'Inferior';
                end
            end
        end
    end
    %%
    % Only make folders for directions we have in the dataset.
    numfov = unique(pixelScale);
    for f=1:length(numfov)
        unique_directions{f} = unique(group_directions(pixelScale==numfov(f)));
    end
    canvas_size = [length(vr) length(ur)];
    
    psconfig( 'pixels', 'pixels', 10, 'no' );
    psnewdoc( canvas_size(2), canvas_size(1), 72, ['SAVE_ME_AS_SOMETHING_NICE.psd'], 'grayscale');
    
    % The sorting goes FOV->Modality->Direction
    for f=length(fovlist):-1:1
        make_Photoshop_group( fovlist{f} );
        for m=MN:-1:1 %Make them backwards so that confocal is on top.
            
            make_Photoshop_group( strrep(ModalitiesSrchStrings{m},'_','') );
            add_to_Photoshop_group( fovlist{f}, 0 );
            
            for k=1: length(unique_directions{f})
                if ~isempty(unique_directions{f}{k})
                    make_Photoshop_group( unique_directions{f}{k} );
                    add_to_Photoshop_group( strrep(ModalitiesSrchStrings{m},'_',''), 0);
                end
                setActiveLayer(fovlist{f}, 1);
            end
        end
        % Make the active layer the next FOV
        setActiveLayer(fovlist{f}, 1);
    end
    
    %%
    for i = 1: NumOfRefs
        for n = RefChains{i}
            loadednames={};
            for m=1: MN
                if ~isempty(imageFilename{m,n})
                    %Read each image, and then transform
                    im = imresize( imread(char(imageFilename{m,n})),pixelScale(n) ); % They need to be pretty for output! Keep it bicubic
                    H = TotalTransform(:,:,n);

                    if size(im,3) == 2
                        im=padarray(im, [length(vr) length(ur) 2]-size(im),0,'post');
                    else
                        im=padarray(im, [length(vr) length(ur)]-size(im),0,'post');
                    end
                    
                    tform = affine2d(H');
                    im_=imwarp(im,tform.invert(),'OutputView',imref2d(size(im)) );
                    
                    [pathstr,name,ext] = fileparts(char(imageFilename{m,n})) ;

                    loadednames{m} = name;

                    if size(im_,3) == 2
                        %add to combined image
                        nonzero = im_(:,:,2)>0;
                        im_ = im_(:,:,1);
                        im_(:,:,1) = uint8(round(im_*255));
                        im_(:,:,2) = uint8(round(nonzero*255));
                    else
                        im_ = im_(:,:,1);                    
                        %add to combined image
                        nonzero = im_>0;
                        im_(:,:,1) = uint8(round(im_*255));
                        im_(:,:,2) = uint8(round(nonzero*255));                         
                    end
                        
                    saveName=[name,'_aligned_to_ref',num2str(i),'_m',num2str(m)];
                    psnewlayer(saveName);

                    pssetpixels(im_(:,:,2),16);
                    pssetpixels(im_(:,:,1), 'undefined');

                    add_to_Photoshop_group(num2str(pixelScale(n),'%0.2f'),1) % FOV
                    add_to_Photoshop_group(ModalitiesSrchStrings{m},0) %Modality
                    add_to_Photoshop_group(group_directions{n},0) % Direction

                    numWritten = numWritten+1;
                    waitbar(numWritten/(N*MN),h,strcat('Writing to Photoshop (',num2str(100*numWritten/(N*MN),3),'%)'));
                end
            end
            % Link the layers we just imported
            link_Photoshop_layers(loadednames);
        end
    end
end

% save tmp.mat;
%%

for m = 1:MN
    
    %initialize blank combined image of all pieces for the modality
    imCombinedAll = zeros(length(vr),length(ur), 'uint8');

    for i = 1: NumOfRefs
        if ~isempty(imageFilename{m,AllRefIndex(i)})
            %initialize blank combined image for the modality/piece
            imCombined = zeros(length(vr),length(ur), 'uint8');

            for n = RefChains{i}
                if ~isempty(imageFilename{m,n})
                    %read each image, and then transform
                    im = imresize( imread(char(imageFilename{m,n})),pixelScale(n) );
                    
                    if size(im,3) == 2
                        im=padarray(im, [size(imCombined) 2]-size(im),0,'post');
                    else
                        im=padarray(im, size(imCombined)-size(im),0,'post');
                    end
                    
                    H = TotalTransform(:,:,n);

                    tform = affine2d(H');
                    im_=imwarp(im,tform.invert(),'OutputView',imref2d(size(im)) );
                    
                    %save each individually transformed image
                    [pathstr,name,ext] = fileparts(char(imageFilename{m,n})) ;

                    if size(im_,3) == 2
                        %add to combined image
                        nonzero = im_(:,:,2)>0;
                        im_ = im_(:,:,1);
                        imCombined(nonzero) = im_(nonzero);
                        im_(:,:,1) = uint8(round(im_*255));
                        im_(:,:,2) = uint8(round(nonzero*255));
                    else
                        im_ = im_(:,:,1);                    
                        %add to combined image
                        nonzero = im_>0;
                        imCombined(nonzero) = im_(nonzero);
                        im_(:,:,1) = uint8(round(im_*255));
                        im_(:,:,2) = uint8(round(nonzero*255));                         
                    end

                    saveFileName=[name,'_aligned_to_ref',num2str(i),'_m',num2str(m),'.tif'];
                    saveTif(im_,outputDir,saveFileName);

                    numWritten = numWritten+1;
                    waitbar(numWritten/(N*MN),h,strcat('Writing Outputs (',num2str(100*numWritten/(N*MN),3),'%)'));

                end
            end
        end
        
        %add to all combined image
        nonzero = imCombined>0;        
        imCombinedAll(nonzero) = imCombined(nonzero);
        imCombined(:,:,2) = round(nonzero*255);
        
        %save combined image for each piece - removed for memory concerns.
        if(NumOfRefs > 1)%only necessary if more than one piece
            saveFileName = strcat('ref_',num2str(i),'_combined_m',num2str(m),'.tif');

            if(AppendToExisting)
                outNameList{i+1,m}=fullfile('Append',saveFileName);
            else
                outNameList{i+1,m}=saveFileName;
            end

            saveTif(imCombined,outputDir,saveFileName);
        end
    end

    %save the combined image of all the pieces
    saveFileName = strcat('all_ref_combined_m',num2str(m),'.tif');
    if(AppendToExisting)
        outNameList{1,m}=fullfile('Append',saveFileName);
    else
        outNameList{1,m}=saveFileName;
    end

    nonzero = imCombinedAll>0;
    imCombinedAll(:,:,2) = round(nonzero*255);

    saveTif(imCombinedAll,outputDir,saveFileName);
    
end

runtime=toc


