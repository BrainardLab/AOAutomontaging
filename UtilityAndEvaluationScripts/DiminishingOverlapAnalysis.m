%Script for evaluating when a pair of overlapping AO images can no longer
%be matched. One column of overlap is removed at each iteration.
%Written by Min Chen (minchen1@upenn.edu)

%Set Input(Raw AO Data) and Output Folders
imageDirBase=getpref('AOAutomontaging','inputOverlapAnalysisDataDir');
outDirBase = getpref('AOAutomontaging','outputAnalysisDir');

imageFileFolders = dir(imageDirBase);
FN=length(imageFileFolders);

 saveAll_new=cell(1,FN);
 complete = zeros(1,FN);

parfor f = find(~complete)%do for every folder
    f
    saveAll_new{f}=cell(5,4,50);
    imageFileFolder          = imageFileFolders(f).name;
    if strcmp(imageFileFolder,'.')||strcmp(imageFileFolder,'..')
        continue;
    end
    
    imageDir=fullfile(imageDirBase,imageFileFolder);
    Allfiles = dir(fullfile(imageDir,'*.tif'));
    Allfiles = {Allfiles.name};
    
    for k = 1:4 %do for each modality option (3 indiv, 1 combined)
        f
        k
        switch k
            case 1
                ModalitiesSrchStrings = {'confocal'};
                MN=1;
                outputDir = 'Confocal';
            case 2
                ModalitiesSrchStrings = {'split_det'};
                MN=1;
                outputDir = 'SplitDetection';
            case 3
                ModalitiesSrchStrings = {'avg'};
                MN=1;
                outputDir = 'DarkField';
            case 4
                ModalitiesSrchStrings = {'confocal'; 'split_det'; 'avg'};
                MN=3;
                outputDir = 'AllThree';
        end
        
        
        
        inData=[];
        
        for m = 1:MN
            if (~isempty(ModalitiesSrchStrings{m}))%check it's not empty
                imagesFound = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, ModalitiesSrchStrings{m}))));%search for image of this modality
                if(~isempty(imagesFound))
                    inData = [inData; imagesFound];
                end
                
            end
        end
        N=2;
        imageFilename = cell(MN,N);%stores all filenames
        for n=1:N
            for m = 1:MN
                imageFilename{m,n} = fullfile(imageDir, inData{m,n});
            end
        end
        
        %stores sift features for each image
        f_all = cell(MN,N);
        d_all = cell(MN,N);
        
        %load image
        im1 = cell(MN,1);
        im2= cell(MN,1);
        for m= 1:MN
            im1{m} = imread(char(imageFilename{m,1}));
            im2{m} = imread(char(imageFilename{m,2}));
        end
        
        stopFlag=0;
        firstItr=1;
        counter=0;
        while(~stopFlag)
            counter=counter+1;
            %Sift features
            for n=1:N
                for m = 1:MN
                    [f1,d1] = vl_sift(im2single(im1{m}),'Levels',55);
                    f_all{m,1} = f1;
                    d_all{m,1} = d1;
                    [f2,d2] = vl_sift(im2single(im2{m}),'Levels',55);
                    f_all{m,2} = f2;
                    d_all{m,2} = d2;
                    
                end
            end
            
            %match
            
            [relativeTransform, numOkMatches, numMatches]=sift_mosaic_fast_MultiModal(im1, im2, [],0,f_all(:,1),d_all(:,1),f_all(:,2),d_all(:,2),1);
            
            %Apply transformation and find overlap(only the first time)
            if( firstItr);
                H=relativeTransform;
                box2 = [1  size(im2{m},2) size(im2{m},2)  1 ;
                    1  1           size(im2{m},1)  size(im2{m},1) ;
                    1  1           1            1 ] ;
                box2_ = inv(H) * box2 ;
                box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
                box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
                ur = min([1 box2_(1,:)]):max([size(im1{m},2) box2_(1,:)]) ;
                vr = min([1 box2_(2,:)]):max([size(im1{m},1) box2_(2,:)]) ;
                
                [u,v] = meshgrid(ur,vr) ;
                im1_ = vl_imwbackward(im2double(im1{m}),u,v) ;
                
                z_ = H(3,1) * u + H(3,2) * v + H(3,3);
                u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
                v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
                im2_ = vl_imwbackward(im2double(im2{m}),u_,v_) ;
                
                mask_=~isnan(im1_) & ~isnan(im2_);
                
                bbox1 = regionprops(~isnan(im1_),'BoundingBox');
                bbox2 = regionprops(~isnan(im2_),'BoundingBox');
                if(bbox1.BoundingBox(1) > bbox1.BoundingBox(2))
                    imToCrop = 1;
                else
                    imToCrop = 2;
                end
                bboxBoth = regionprops(mask_,'BoundingBox');
                pixelColLeft=bboxBoth.BoundingBox(3);
                colscropped=0;
                firstItr=0;
            end
            %Record transformation, numOkMatches, numMatches, number of pixels columns
            %overlapping, @ modality, folder number
            
            saveAll_new{f}{1,k,counter}=relativeTransform;
            saveAll_new{f}{2,k,counter}=numOkMatches;
            saveAll_new{f}{3,k,counter}=numMatches;
            saveAll_new{f}{4,k,counter}=pixelColLeft;
            saveAll_new{f}{5,k,counter}=colscropped;
            
            %crop Image 2, unless it's less than a column, in which case we stop.
            if(pixelColLeft > 400)
                colscropped=colscropped+pixelColLeft-400;
                pixelColLeft=400;
            elseif(pixelColLeft > 300)
                colscropped=colscropped+pixelColLeft-300;
                pixelColLeft=300;
            elseif(mod(pixelColLeft,10) ~= 0)
                colscropped=colscropped+mod(pixelColLeft,10);
                pixelColLeft=pixelColLeft-mod(pixelColLeft,10);
            elseif(pixelColLeft > 200)
                colscropped=colscropped+10;
                pixelColLeft=pixelColLeft-10;
            elseif(mod(pixelColLeft,5) ~= 0)
                colscropped=colscropped+mod(pixelColLeft,5);
                pixelColLeft=pixelColLeft-mod(pixelColLeft,5);
            else
                colscropped=colscropped+5;
                pixelColLeft=pixelColLeft-5;
            end
            
            
            if(imToCrop == 1)
                
                for m = 1:MN
                    temp=im1{m};
                    temp(:,1:colscropped)=0;
                    im1{m}=temp;
                end
            else
                for m = 1:MN
                    temp=im2{m};
                    temp(:,1:colscropped)=0;
                    im2{m}=temp;
                end
            end
            
            if(pixelColLeft<=0 || numOkMatches <=2)
                stopFlag=1;
            end
        end
    end
    %repeat
    complete(f) = 1;
end

%save results
save('OverlapAnalysis.mat','saveAll','-v7.3');
