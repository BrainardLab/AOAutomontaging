function [bestAlignTransMtx, numOkMatches_all, numMatches_all, bestScale]= matchImagesMultiModalFeatures(im1, im2, saveDir,saveFlag,f1,d1,f2,d2,TransType,rotLimit,featureType,modalitiesToUse)
% Matches two images using given precalculated SIFT features
%Input:
%im1 -- Reference Image
%im2 -- Moving Image to be matched
%saveMatchesName -- String name if saving the matches to a figure
%saveFlag -- Set to 1 if you would like to generate figure showing matches
%f1, f2 -- Vectors of SIFT feature locations locations for im1 and im2
%d1, d2 -- Vectors of SIFT feature descriptors corresponding to f1 and f2
%TransType -- Index for the type of transformation used by the matching:
%  0 - Translation Only
%  1 - Translation + Rotation
%  2 - Rotation + Translation + Scaling
%  3 - Affine
%  4 - Rotation + Translation + Individual Scaling
%rotLimit -- The maximum rotation allowed in the transformation
%featureType -- Feature type to use for matching
%  0 - SIFT
%  1 - Constellation Feature
%  2 - SURF

%Outputs:
%bestAlignTransMtx -- Best transformation found (in matrix form)
%numOkMatches_all -- Total number of inlier matches determined by RANSAC
%numMatches_all -- Total number of matches found
%BestScale -- Best scale found for the transform

%Written by Min Chen (minchen1@upenn.edu)

% --------------------------------------------------------------------
%                                                         SIFT matches
% --------------------------------------------------------------------

%default limit for roation is 10 degrees
if ~exist('rotLimit','var') || isempty(rotLimit)
    rotLimit=pi/18;
end

%default feature type is SIFT
if ~exist('featureType','var') || isempty(featureType)
    featureType=0;
end

%parameter
matchTolerance = 6;%determines how close a tansformed match needs to be to be considerd an inlier
if(featureType==0 || featureType==2)%change this depending on the type of feature
    matchTolerance = 6;%sift has some leeway due to filtering
elseif(featureType==1)
    matchTolerance = 6;%10;%Constellation features should be close to cell-to-cell
end
%saveFlag=1;
%saveDir='C:\Users\dontm\Downloads\temp\New folder (7)';

MN = size(f1,1);
MRange=1:MN;%modalities to include in RANSAC

if(exist('modalitiesToUse','var') && find(modalitiesToUse,1,'last')<=MN)
    MRange = find(modalitiesToUse,1,'first'):find(modalitiesToUse,1,'last');
end
%check if either image had zero features, exit with no matches if so
total_f1 = 0;
total_f2 = 0;
for m = 1:MN
    total_f1 = total_f1+size(f1{m},2);
    total_f2 = total_f2+size(f2{m},2);
end
if(total_f1 < 1 || total_f2 < 1)
    bestAlignTransMtx = eye(3,3);
    numOkMatches_all=0;
    numMatches_all=0;
    bestScale=1;
    return
end



X1 = [];
X2 = [];
matches = cell(MN,1);
numMatches = zeros(MN,1);
for m = MRange
    
    %use differen matching method depending on feature type
    if(featureType==0)%sift
        matches_m = matchFeatures(d1{m}',d2{m}',MatchThreshold=100);
        matches_m = matches_m';%transpose for consistency with rest of code
        %[matches_m, scores] = vl_ubcmatch_fast(d1{m},d2{m});
    elseif(featureType==1)%constellation
        matches_m = matchGridFeatures(d1{m},d2{m});%,d1_int,d2_int);
    elseif(featureType==2)%surf
        matches_m = matchFeatures(d1{m},d2{m});
        matches_m = matches_m';%transpose for consistency with rest of code
        f1{m} = f1{m}.Location';
        f2{m} = f2{m}.Location';
    end
    X1_m = f1{m}(1:2,matches_m(1,:)) ;
    X2_m = f2{m}(1:2,matches_m(2,:)) ;
    
    
    %check for duplicate matches
    %check duplicates in X1 & X2 pairs
    [uX1_m, IA, IC1] = unique([round(X1_m') round(X2_m')],'rows');
    X1_m = X1_m(1:2,IA');
    X2_m = X2_m(1:2,IA');
    matches_m = matches_m(:,IA');
    
    matches{m} = matches_m;
    numMatches(m) = size(X1_m,2);
    X1=[X1 X1_m];
    X2=[X2 X2_m];
end

%Homogeneous coordinates of a 2D point, so 3rd element is 1
X1(3,:) = 1 ;
X2(3,:) = 1 ;


%numMatches = size(matches,2) ;
numMatches_all = size(X1,2);

% --------------------------------------------------------------------
%                                         RANSAC
% --------------------------------------------------------------------

clear H score ok ;
bestScore = -1;
bestAlignTransMtx = eye(3,3);
bestOK_all = zeros(1,size(X1,2));
bestScale = 1;

allIndex = 1:numMatches_all;%list of all the indexs for the matches
offset = 0;
for m = [0 MRange]
    
    if (m==0)
        testIndices = allIndex;%the first test is across all modalities
    else%then we start testing individual modalities
        testRange = zeros(1,length(allIndex));
        testRange(offset+1:offset+numMatches(m))=1;
        offset = offset + numMatches(m);
        testIndices = allIndex(logical(testRange));
    end
    for t = 1:10000
        
        % estimate model
        if(length(testIndices) >= 3)
            subset = randsample(testIndices,3);
            %subset = vl_colsubset(1:numMatches_all, 3) ;
        else
            subset = testIndices;
        end
        
        H = eye(3,3);
        TransRadians = 0;
        Scale = 1;
        if (TransType == 0)
            %Score Using Only Translation
            C1 = mean(X1(1:2,subset),2);
            C2 = mean(X2(1:2,subset),2);
            trans = C2 - C1;
            H(1,3) = trans(1);
            H(2,3) = trans(2);
            TransRadians = 0;
        elseif (TransType == 1)
            %Score Using Only rotation + translation
            [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', false, 'reflection', false);
            H = [trans.T' mean(trans.c,1)'; 0 0 1];
            TransRadians = atan2(H(1,2),H(1,1));
        elseif (TransType == 2)
            %Score Using rotation + translation + scaling
            [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', true, 'reflection', false);
            H = [trans.b*trans.T' mean(trans.c,1)'; 0 0 1];
            T = trans.T';
            TransRadians = atan2(T(1,2),T(1,1));
            Scale = trans.b;
        elseif(TransType == 3)
            %Score using homography (NOT affine)
            %         A = [] ;
            %         for i = subset
            %             A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
            %         end
            %         [U,S,V] = svd(A) ;
            %         H = reshape(V(:,9),3,3) ;
            
            % Use MATLAB function until redoing with proper affine (parallel
            % lines must be preserved!)
            trans = estimateGeometricTransform(X1(1:2,:)', X2(1:2,:)','affine', 'MaxNumTrials',1000,'MaxDistance',matchTolerance);
            H = trans.T';
            
            TransRadians = atan2(H(1,2),H(1,1));
        elseif(TransType == 4)
            %individual scale
            [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', false, 'reflection', false);
            H = [trans.T' mean(trans.c)'; 0 0 1];
            a=X2(1,subset)*pinv(Z(:,1)');
            b=X2(2,subset)*pinv(Z(:,2)');
            H = [a 0 0; 0 b 0; 0 0 1] * H;
            TransRadians = atan2(H(1,2),H(1,1));
        end
        
        
        X2_ = H * X1 ;
        du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
        dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
        ok = (du.*du + dv.*dv) < matchTolerance*matchTolerance;
        score = sum(ok);
        
        if((score > bestScore) && ((abs(TransRadians) < rotLimit)) && Scale < 1.1 && Scale > .90)
            bestScore = score;
            bestAlignTransMtx = H;
            bestOK_all = ok;
            bestScale = Scale;
        end
        
    end
end
%bestScale;
numOkMatches_all = sum(bestOK_all);
numOkMatches = zeros(MN,1);
%separate valid matches by modality (visualization purpose only)
bestOK = cell(MN,1);
offset = 0;
for m = 1:MN
    bestOK{m} = bestOK_all(offset+1:offset+numMatches(m));
    numOkMatches(m) = sum(bestOK{m});
    offset = offset + numMatches(m);
end

% --------------------------------------------------------------------
%                                                Show and save matches
% --------------------------------------------------------------------
if(saveFlag>=1)
    colorRGB = [51/255 153/255 255/255];
    colorRGB2 = [255/255 0/255 0/255];
    colorRGB3 = [0/255 255/255 0/255];
    
    figID=figure(1); clf;
    title(sprintf('%d total tentative matches', numMatches_all)) ;
    title(sprintf('%d (%.2f%%) total inliner matches out of %d', ...
        sum(bestOK_all), ...
        100*sum(bestOK_all)/numMatches_all, ...
        numMatches_all)) ;
    
    %do for every modality
    for m = MRange
        dh1 = max(size(im2{m},1)-size(im1{m},1),0) ;
        dh2 = max(size(im1{m},1)-size(im2{m},1),0) ;

        %plot all matches
        figID1 = figure(1);
        gapSize = 25;
        imOut = [padarray(im1{m}(:,:,1),[dh1 gapSize],255,'post') padarray(im2{m}(:,:,1),dh2,255,'post')];
        imagesc(imOut);
        caxis([min(min(min(im1{m}(:,:,1))),min(min(im2{m}(:,:,1)))) max(max(max(im1{m}(:,:,1))),max(max(im2{m}(:,:,1))))]);
        o = size(im1{m},2) +gapSize;
        if(numMatches(m) > 0)
            line([f1{m}(1,matches{m}(1,:));f2{m}(1,matches{m}(2,:))+o], ...
                [f1{m}(2,matches{m}(1,:));f2{m}(2,matches{m}(2,:))],'color',colorRGB,'LineWidth',2) ;
        end
        title(sprintf('%d tentative matches', numMatches(m))) ;
        axis image off ;
        colormap('gray')
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveMatchesName=['AllMatches_m' num2str(m) '.bmp'];
        saveas(figID1,fullfile(saveDir,saveMatchesName))
        
        %plot inlier matches
        figID2 = figure(1);
        imOut = [padarray(im1{m}(:,:,1),[dh1 gapSize],255,'post') padarray(im2{m}(:,:,1),dh2,255,'post')];
        imagesc(imOut);
        caxis([min(min(min(im1{m}(:,:,1))),min(min(im2{m}(:,:,1)))) max(max(max(im1{m}(:,:,1))),max(max(im2{m}(:,:,1))))]);
        o = size(im1{m},2) + gapSize;
        if(max(bestOK{m}) > 0)
            line([f1{m}(1,matches{m}(1,bestOK{m}));f2{m}(1,matches{m}(2,bestOK{m}))+o], ...
                [f1{m}(2,matches{m}(1,bestOK{m}));f2{m}(2,matches{m}(2,bestOK{m}))],'color',colorRGB,'LineWidth',2) ;
        end
        title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
            sum(bestOK{m}), ...
            100*sum(bestOK{m})/numMatches(m), ...
            numMatches(m))) ;
        axis image off ;
        colormap('gray')
        
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveMatchesName=['RANSACMatches_m' num2str(m) '.bmp'];
        saveas(figID2,fullfile(saveDir,saveMatchesName))
        
        %plot location of features
        figID2 = figure(1);
        imOut = [padarray(im1{m}(:,:,1),[dh1 gapSize],255,'post') padarray(im2{m}(:,:,1),dh2,255,'post')];
        imagesc(imOut);
        caxis([min(min(min(im1{m}(:,:,1))),min(min(im2{m}(:,:,1)))) max(max(max(im1{m}(:,:,1))),max(max(im2{m}(:,:,1))))]);
        o = size(im1{m},2) + gapSize;
        hold on
        scatter([f1{m}(1,:) (f2{m}(1,:)+o)], ...
            [f1{m}(2,:) (f2{m}(2,:))],2,colorRGB3) ;
        
%         if(max(bestOK{m}) > 0)
%             scatter([f1{m}(1,matches{m}(1,bestOK{m})) (f2{m}(1,matches{m}(2,bestOK{m}))+o)], ...
%                 [f1{m}(2,matches{m}(1,bestOK{m})) f2{m}(2,matches{m}(2,bestOK{m}))],2,colorRGB2) ;
%         end
        
        hold off
        title(sprintf('feature locations')) ;
        axis image off ;
        colormap('gray')
        
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveMatchesName=['FeatureLocations_m' num2str(m) '.bmp'];
        saveas(figID2,fullfile(saveDir,saveMatchesName))
    end
    
    
end

% --------------------------------------------------------------------
%                                                 Show Matched  Mosaic
% --------------------------------------------------------------------
%
if(saveFlag>=1)
    H = bestAlignTransMtx;
    %do for 
    for m = MRange
        box2 = [1  size(im2{m},2) size(im2{m},2)  1              ;
            1  1              size(im2{m},1)  size(im2{m},1) ;
            1  1              1               1            ] ;
        box2_ = inv(H) * box2 ;
        box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
        box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
        ur = min([1 box2_(1,:)]):max([size(im1{m},2) box2_(1,:)]) ;
        vr = min([1 box2_(2,:)]):max([size(im1{m},1) box2_(2,:)]) ;
        
        [u,v] = meshgrid(ur,vr) ;
        im1_ = interp2(im2double(im1{m}(:,:,1)),u,v);
        im1_ = im1_(:,:,1);
        
        saveTif(im1_,saveDir,['image_A_m' num2str(m) '.tif']);
        
        z_ = H(3,1) * u + H(3,2) * v + H(3,3);
        u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
        v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
        im2_ = interp2(im2double(im2{m}(:,:,1)),u_,v_);
        im2_ = im2_(:,:,1);
        
        saveTif(im2_,saveDir,['image_B_Trans_m' num2str(m) '.tif']);
        
        if(saveFlag==1)
            im1_(isnan(im1_)) = 0 ;
            im2_(isnan(im2_)) = 0 ;
            
            figID4 = figure(4);
            imshowpair(im1_,im2_) ; axis image off ;
            set(gca,'LooseInset',get(gca,'TightInset'));
            saveMatchesName=['compare_m' num2str(m) '.bmp'];
            saveas(figID4,fullfile(saveDir,saveMatchesName))
                        
            overlap= (im1_ ~= 0) & (im2_ ~= 0);
            imAvg = (im1_+im2_)./(overlap+ones(size(overlap)));
            saveTif(imAvg,saveDir,['image_Avg_m' num2str(m) '.tif']);
            
            im1_temp = im1_;
            im1_temp(overlap) = 0; %for visualization, im2 is ontop
            mosaic = (im1_temp + im2_);
            
            figID3 = figure(3) ; clf ;
            imagesc(mosaic) ; axis image off ;
            colormap('gray')
            
            set(gca,'LooseInset',get(gca,'TightInset'));
            saveMatchesName=['mosaic_m' num2str(m) '.bmp'];
            saveas(figID3,fullfile(saveDir,saveMatchesName))
            
            hold on
            f2_ = inv(H)*[f2{m}(1,:); f2{m}(2,:); ones(size(f2{m}(1,:)))];
            scatter(f1{m}(1,:)-ur(1)+1, f1{m}(2,:)-vr(1)+1,2,colorRGB) ;
            scatter(f2_(1,:)-ur(1)+1, f2_(2,:)-vr(1)+1,2,colorRGB3) ;
            
            if(max(bestOK{m}) > 0)
                f2OK_ = inv(H)*[f2{m}(1,matches{m}(2,bestOK{m})); f2{m}(2,matches{m}(2,bestOK{m})); ones(size(f2{m}(2,matches{m}(2,bestOK{m}))))];
                scatter(f1{m}(1,matches{m}(1,bestOK{m}))-ur(1)+1, ...
                    f1{m}(2,matches{m}(1,bestOK{m}))-vr(1)+1,2,colorRGB2) ;
                scatter(f2OK_(1,:)-ur(1)+1, f2OK_(2,:)-vr(1)+1,2,colorRGB2) ;
            end
            hold off
            
            saveMatchesName=['mosaicFeatures_m' num2str(m) '.bmp'];
            saveas(figID3,fullfile(saveDir,saveMatchesName))
            
            %save individual with features
            figID3 = figure(3) ; clf ;
            imagesc(im1_) ; axis image off ;
            colormap('gray')
            
            set(gca,'LooseInset',get(gca,'TightInset'));
            hold on
            scatter(f1{m}(1,:)-ur(1)+1, f1{m}(2,:)-vr(1)+1,2,colorRGB) ;
            
            if(max(bestOK{m}) > 0)
                scatter(f1{m}(1,matches{m}(1,bestOK{m}))-ur(1)+1, ...
                    f1{m}(2,matches{m}(1,bestOK{m}))-vr(1)+1,2,colorRGB2) ;
            end
            hold off
            
            saveMatchesName=['image_Afeatures_m' num2str(m) '.bmp'];
            saveas(figID3,fullfile(saveDir,saveMatchesName))
            
            figID3 = figure(3) ; clf ;
            imagesc(im2_) ; axis image off ;
            colormap('gray')
            
            set(gca,'LooseInset',get(gca,'TightInset'));
            hold on
            scatter(f2_(1,:)-ur(1)+1, f2_(2,:)-vr(1)+1,2,colorRGB3) ;
            
            if(max(bestOK{m}) > 0)
                scatter(f2OK_(1,:)-ur(1)+1, f2OK_(2,:)-vr(1)+1,2,colorRGB2) ;
            end
            hold off
            
            saveMatchesName=['image_B_Transfeatures_m' num2str(m) '.bmp'];
            saveas(figID3,fullfile(saveDir,saveMatchesName))
        end
    end
end

