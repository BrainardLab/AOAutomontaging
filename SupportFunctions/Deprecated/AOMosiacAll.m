function  outNameList = AOMosiacAll(imageDir,posFileLoc,outputDir)

%Load Data
%baseDir = 'C:\Users\Min\Dropbox\AOSLOImageProcessing\ConstructMontage\InputImages_Set1\';
%fileID = fopen(fullfile(baseDir,'ImageInfo.txt'));

%Load info from descriptor file
% formatSpec = '%s';
% C_text = textscan(fileID,formatSpec,9,'Delimiter','\t');
% C = textscan(fileID, '%d8 %s %s %s %s %d8 %s %d8 %d8', 'delimiter', '\t');
%N = size(C{1},1);

%imageDir = 'C:\Users\Min\Dropbox\AOSLOImageProcessing\ConstructMontage\InputImages_Set3\Images';
%Allfiles = dir('C:\Users\Min\Dropbox\AOSLOImageProcessing\ConstructMontage\InputImages_Set3\Images\*.tif');
Allfiles = dir(fullfile(imageDir,'*.tif'));
Allfiles = {Allfiles.name};
confocal = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_confocal_'))));
darkfield = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_avg_'))));
splitDecision = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_split_'))));

inData = confocal;

N = size(inData,2);

%initialization
LocXY = 100000*ones(2,N);%set as a ridiculous coordinate, incase can't find in spreadsheet
RelativeTransformToRef = zeros(3,3,N);
imageFilename = cell(N,1);
Matched = zeros(1,N);
MatchedTo = zeros(1,N);

ResultsNumOkMatches = -ones(N,N);
ResultsNumMatches = -ones(N,N);
ResultsTransformToRef = zeros(3,3,N,N);

%load info from excel
%Dataset 2
%[temp,C,temp] = xlsread('C:\Users\Min\Dropbox\AOSLOImageProcessing\ConstructMontage\InputImages_Set2\ImageInfoRaw_CS_13213_20150302_OS.xlsx','A9:C87');
%Dataset 3
%'C:\Users\Min\Documents\Research\AdaptiveOpticsMosaic\NewDataSet3_2015_8_25\ImageInfoRaw_NC_11028_20140527.xlsx'
[temp,C,temp] = xlsread(posFileLoc,'A1:C59');

for n=1:N
    imageFilename{n} = fullfile(imageDir, inData{n});
    for i = 1:size(C,1)
        if (~isempty(strfind(inData{n}, C{i,1})))
            Loc = strsplit(C{i,3},',');
            LocXY(1,n) = str2double(strtrim(Loc{1}));
            LocXY(2,n) = str2double(strtrim(Loc{2}));
            break;
        end
    end
end
   

%old format
% for n=1:N
%    imageFilename{n} = fullfile(baseDir, strcat(C{2}(n),C{5}(n),C{3}(n),num2str(C{6}(n)),C{4}(n),'.tif'));
%    LocXY(1,n) = C{8}(n);
%    LocXY(2,n) = C{9}(n);
% end
%[indexStructs(i).prefix indexStructs(i).videoAndRef indexStructs(i).middle num2str(indexStructs(i).nNum) indexStructs(i).suffix '.tif']);

%sort using LocXY
[sortLocXY I] = sortrows(abs(LocXY)',[1,2]);
LocXY=LocXY(:,I);
imageFilename=imageFilename(I);

%Calculate All Sift Features
f_all = cell(N,1);
d_all = cell(N,1);

h = waitbar(0,'Calculating Sift Features 0%');
for n=1:N
im = im2single(imread(char(imageFilename{n})));
[f1,d1] = vl_sift(im,'Levels',55);
f_all{n} = f1;
d_all{n} = d1;
waitbar(n/N,h,strcat('Calculating Sift Features (',num2str(100*n/N,3),'%)'));
end

%Begin Matching
while (sum(Matched) < N)

    waitbar(sum(Matched)/N,h,strcat('Aligning Images (',num2str(100*sum(Matched)/N,3),'%)'));
    
    %find the closest unmatched image to the origin. Set as new reference frame
    unMatchedI = find(Matched == 0);
    unMatchedLocXY = LocXY(:,unMatchedI);
    [sortUnMatchedLocXY I] = sortrows(abs(unMatchedLocXY)',[1,2]);
    
    %set reference frame
    refIndex = unMatchedI(I(1))
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
        for n=2:N
            if(Matched(n) == 0) % look for an unmatached Images
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

                            if (ResultsNumOkMatches(n,refIndex) == -1) %skip if match has already been found before

                                refImg = imread(char(imageFilename{refIndex}));
                                currentImg = imread(char(imageFilename{n}));

                                saveMatchesName = strcat('Matches_n',num2str(n),'_(', num2str(LocXY(1,n)),',', ...
                                num2str(LocXY(2,n)),')','_to_','n',num2str(refIndex),'_(', num2str(LocXY(1,refIndex)),...
                                ',',num2str(LocXY(2,refIndex)),')','.jpg');
                                [relativeTransform, numOkMatches, numMatches]=sift_mosaic_fast(refImg, currentImg, saveMatchesName,1,f_all{refIndex},d_all{refIndex},f_all{n},d_all{n},0);

                                ResultsNumOkMatches(n,refIndex) = numOkMatches;
                                ResultsNumMatches(n,refIndex) = numMatches;
                                ResultsTransformToRef(:,:,n,refIndex) = relativeTransform;
                            end

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

outNameList = cell(NumOfRefs,1);

for i = 1: NumOfRefs
    %find max and min bounding box over all iamges attached to this ref

    maxX =-1000000000; 
    minX =1000000000;
    maxY =-1000000000; 
    minY =1000000000;
    
    for n = 1:N
        im = imread(char(imageFilename{n}));
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
    
    im =  imread(char(imageFilename{AllRefIndex(i)}));
    imCombined = vl_imwbackward(im2double(im),u,v);

    for n = RefChains{i}
        im = imread(char(imageFilename{n}));
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
        [pathstr,name,ext] = fileparts(char(imageFilename{n})) ;
        imwrite(im_,fullfile(outputDir,[name,'_aligned.png']),'alpha',double(~isnan(im_))); 
    end

    %figure(i);clf;
    %imshow(imCombined)
    outName = strcat('ref_',num2str(i),'_combined','.png');
    outNameList{i}=outName;
    imwrite(imCombined,fullfile(outputDir,outName),'alpha',double(~isnan(imCombined))); 
end