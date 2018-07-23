function [bestH, numOkMatches_all, numMatches_all, bestScale]= sift_mosaic_fast_MultiModal(im1, im2, saveDir,saveFlag,f1,d1,f2,d2,TransType,rotLimit,featureType)
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

%Outputs:
%bestH -- Best transformation found (in matrix form)
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

if ~exist('featureType','var') || isempty(featureType)
    featureType=0;
end

%parameter
matchTolerance = 6;%determines how close a tansformed match needs to be to be considerd an inlier
if(featureType==0)%change this depending on the type of feature
    matchTolerance = 6;%sift has some leeway due to filtering
elseif(featureType==1)
    matchTolerance = 3;%Constellation features should be close to cell-to-cell
end
%saveFlag=1;
%saveDir='C:\Users\dontm\Downloads\temp\New folder (7)';

MN = size(f1,1);

%check if either image had zero features, exit with no matches if so
total_f1 = 0;
total_f2 = 0;
for m = 1:MN
    total_f1 = total_f1+size(f1{m},2);
    total_f2 = total_f2+size(f2{m},2);
end
if(total_f1 < 1 || total_f2 < 1)
    bestH = eye(3,3);
    numOkMatches_all=0;
    numMatches_all=0;
    bestScale=1;
    return
end



X1 = [];
X2 = [];
matches = cell(MN,1);
numMatches = zeros(MN,1);
for m = 1:MN
    
    if(featureType==0)
        [matches_m, scores] = vl_ubcmatch_fast(d1{m},d2{m});
    elseif(featureType==1)
        matches_m = matchGridFeatures(d1{m},d2{m});%,d1_int,d2_int);
    end
    
    X1_m = f1{m}(1:2,matches_m(1,:)) ;
    X2_m = f2{m}(1:2,matches_m(2,:)) ;
    %check for duplicate matches
    %check duplicates in X1 & X2 pairs
    [uX1_m, IA, IC1] = unique([round(X1_m') round(X2_m')],'rows');
    X1_m = X1_m(1:2,IA');
    X2_m = X2_m(1:2,IA');
    matches_m = matches_m(:,IA');
    
    %check X2 next
    %    [uX2_m, IA, IC1] = unique(round(X2_m'),'rows');
    %    X1_m = X1_m(1:2,IA');
    %    X2_m = X2_m(1:2,IA');
    %    matches_m = matches_m(:,IA');
    
    matches{m} = matches_m;
    numMatches(m) = size(X1_m,2);
    X1=[X1 X1_m];
    X2=[X2 X2_m];
end


X1(3,:) = 1 ;
X2(3,:) = 1 ;


%numMatches = size(matches,2) ;
numMatches_all = size(X1,2);

% --------------------------------------------------------------------
%                                         RANSAC with homography model
% --------------------------------------------------------------------

clear H score ok ;
bestScore = -1;
bestH = eye(3,3);
bestOK_all = zeros(1,size(X1,2));
bestScale = 1;

allIndex = 1:numMatches_all;%list of all the indexs for the matches
for t = 1:5000
    % estimate model
    if(numMatches_all >= 3)
        subset = randsample(allIndex,3);
        %subset = vl_colsubset(1:numMatches_all, 3) ;
    else
        subset = allIndex;
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
    score = sum(ok) ;
    
    if((score > bestScore) && ((abs(TransRadians) < rotLimit)) && Scale < 1.1 && Scale > .90)
        bestScore = score;
        bestH = H;
        bestOK_all = ok;
        bestScale = Scale;
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
%                                                  Optional refinement
% --------------------------------------------------------------------

% function err = residual(H)
%  u = H(1) * X1(1,ok) + H(4) * X1(2,ok) + H(7) ;
%  v = H(2) * X1(1,ok) + H(5) * X1(2,ok) + H(8) ;
%  d = H(3) * X1(1,ok) + H(6) * X1(2,ok) + 1 ;
%  du = X2(1,ok) - u ./ d ;
%  dv = X2(2,ok) - v ./ d ;
%  err = sum(du.*du + dv.*dv) ;
% end
%
% if exist('fminsearch') == 2
%   H = H / H(3,3) ;
%   opts = optimset('Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8) ;
%   H(1:8) = fminsearch(@residual, H(1:8)', opts) ;
% else
%   warning('Refinement disabled as fminsearch was not found.') ;
% end

% --------------------------------------------------------------------
%                                                         Show matches
% --------------------------------------------------------------------
if(saveFlag)
    %[saveDir,name,ext] = fileparts(saveMatchesName);
     colorRGB = [51/255 153/255 255/255];
     f1_ok = f1{1}(1:2,matches{1}(1,bestOK{1}));
     f2_ok = f2{1}(1:2,matches{1}(2,bestOK{1}));
     d1_ok = d1{1}(1:2,matches{1}(1,bestOK{1}));
     d2_ok = d2{1}(1:2,matches{1}(2,bestOK{1}));
     colorRGB2 = [255/255 153/255 0/255];
     k=9
    figID=figure(1) ; clf ;
    title(sprintf('%d total tentative matches', numMatches_all)) ;
    title(sprintf('%d (%.2f%%) total inliner matches out of %d', ...
        sum(bestOK_all), ...
        100*sum(bestOK_all)/numMatches_all, ...
        numMatches_all)) ;
    
    
    for m = 1:MN
        dh1 = max(size(im2{m},1)-size(im1{m},1),0) ;
        dh2 = max(size(im1{m},1)-size(im2{m},1),0) ;
        %plot all matches
        %subplot(2,MN,m) ;
        figID1 = figure(1);
        gapSize = 25;
        imOut = [padarray(im1{m}(:,:,1),[dh1 gapSize],255,'post') padarray(im2{m}(:,:,1),dh2,255,'post')];
        imagesc(imOut);
        caxis([min(min(min(im1{m}(:,:,1))),min(min(im2{m}(:,:,1)))) max(max(max(im1{m}(:,:,1))),max(max(im2{m}(:,:,1))))]);
        o = size(im1{m},2) +gapSize;
        if(numMatches(m) > 0)
            line([f1{m}(1,matches{m}(1,:));f2{m}(1,matches{m}(2,:))+o], ...
                [f1{m}(2,matches{m}(1,:));f2{m}(2,matches{m}(2,:))],'color',colorRGB,'LineWidth',2) ;
           % line([f1_ok(1,k);f2_ok(1,k)+o], ...
            %    [f1_ok(2,k);f2_ok(2,k)],'color',colorRGB2,'LineWidth',2) ;

        end
        title(sprintf('%d tentative matches', numMatches(m))) ;
        axis image off ;
        colormap('gray')
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveMatchesName=['SIFTmatches_m' num2str(m) '.bmp'];
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
           % line([f1_ok(1,k);f2_ok(1,k)+o], ...
            %    [f1_ok(2,k);f2_ok(2,k)],'color',colorRGB2,'LineWidth',2) ;

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
        
        %plot select areas
%         figID3 = figure(1);
%         imOut = [padarray(im1{m}(:,:,1),[dh1 gapSize],255,'post') padarray(im2{m}(:,:,1),dh2,255,'post')];
%         imagesc(imOut);
%         caxis([min(min(min(im1{m})),min(min(im2{m}))) max(max(max(im1{m})),max(max(im2{m})))]);
%         o = size(im1{m},2) + gapSize;
%         w = 75;
%         h = 50;
%         if(max(bestOK{2}) > 0)
%             
%             f1_ok = f1{2}(1:2,matches{2}(1,bestOK{2}));
%             f2_ok = f2{2}(1:2,matches{2}(2,bestOK{2}));
%             
%             for q = [3 10] 1:length(f1_ok)
%                 q
%                 if(q == 3)
%                     colorRGB = [255/255 153/255 0/255];
%              %   elseif(q == 12)
%              %       colorRGB = [255/255 0/255 0/255];
%                 else
%                   colorRGB = [51/255 153/255 255/255];
%                 end
%                 posA = [f1_ok(1,q)-w f1_ok(2,q)-h 2*w+1 2*h+1];
%                 rectangle('Position',posA,'EdgeColor',colorRGB,'LineWidth',2);
%                 x_range = round(f1_ok(1,q)-w):round(f1_ok(1,q)+w);
%                 y_range = round(f1_ok(2,q)-h):round(f1_ok(2,q)+h);
%                 im_recA = imOut(y_range,x_range);
%                 posB = [f2_ok(1,q)+o-w f2_ok(2,q)-h 2*w+1 2*h+1];
%                 rectangle('Position',posB,'EdgeColor',colorRGB,'LineWidth',2);
%                 x_range = round(f2_ok(1,q)-w+o):round(f2_ok(1,q)+w+o);
%                 y_range = round(f2_ok(2,q)-h):round(f2_ok(2,q)+h);
%                 im_recB = imOut(y_range,x_range);
%                 imwrite(im_recA,fullfile(saveDir,['imA_rect' num2str(q) '_m' num2str(m) '_.tif']))
%                 imwrite(im_recB,fullfile(saveDir,['imB_rect' num2str(q) '_m' num2str(m) '_.tif']))
% 
%                 if (m == 4)
%                     imA_grid=convertToGrid(im_recA,w,5);
%                     imB_grid=convertToGrid(im_recB,w,5);
%                     imwrite(imA_grid,fullfile(saveDir,['imA_rect' num2str(q) '_m' num2str(m) '_grid_.tif']))
%                     imwrite(imB_grid,fullfile(saveDir,['imB_rect' num2str(q) '_m' num2str(m) '_grid_.tif']))
%                     figure(5);clf
%                     imshowpair(imA_grid,imB_grid);
%                 end                
%                 %            line([f1{m}(1,matches{m}(1,bestOK{m}));f2{m}(1,matches{m}(2,bestOK{m}))+o], ...
%                 %                [f1{m}(2,matches{m}(1,bestOK{m}));f2{m}(2,matches{m}(2,bestOK{m}))]) ;
%             end
%         end
%         title('matching areas');
%         axis image off ;
%         colormap('gray')
%         
%         set(gca,'LooseInset',get(gca,'TightInset'));
%         saveMatchesName=['MatchedAreas_m' num2str(m) '.bmp'];
%         saveas(figID2,fullfile(saveDir,saveMatchesName))
        
     end
    
    
end

% --------------------------------------------------------------------
%                                                        Show   Mosaic
% --------------------------------------------------------------------
%
if(saveFlag)
    H = bestH;
    for m = 1:MN
        box2 = [1  size(im2{m},2) size(im2{m},2)  1              ;
            1  1              size(im2{m},1)  size(im2{m},1) ;
            1  1              1               1            ] ;
        box2_ = inv(H) * box2 ;
        box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
        box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
        ur = min([1 box2_(1,:)]):max([size(im1{m},2) box2_(1,:)]) ;
        vr = min([1 box2_(2,:)]):max([size(im1{m},1) box2_(2,:)]) ;
        
        [u,v] = meshgrid(ur,vr) ;
        im1_ = vl_imwbackward(im2double(im1{m}(:,:,1)),u,v) ;
        im1_ = im1_(:,:,1);
        
        saveTif(im1_,saveDir,['image_A_m' num2str(m) '.tif']);
        
        z_ = H(3,1) * u + H(3,2) * v + H(3,3);
        u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
        v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
        im2_ = vl_imwbackward(im2double(im2{m}(:,:,1)),u_,v_) ;
        im2_ = im2_(:,:,1);
        
        saveTif(im2_,saveDir,['image_B_Trans_m' num2str(m) '.tif']);
        
        
        % mass = ~isnan(im1_) + ~isnan(im2_) ;
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
        
        
        im1_(overlap) = 0; %for visualization, im2 is ontop
        mosaic = (im1_ + im2_);
        
        figID3 = figure(3) ; clf ;
        imagesc(mosaic) ; axis image off ;
        colormap('gray')
        %title('Mosaic') ;
        
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveMatchesName=['mosaic_m' num2str(m) '.bmp'];
        saveas(figID3,fullfile(saveDir,saveMatchesName))
        
    end
end

