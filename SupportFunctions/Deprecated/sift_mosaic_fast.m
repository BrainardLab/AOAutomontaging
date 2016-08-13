function [bestH, numOkMatches, numMatches]= sift_mosaic_fast(im1, im2, saveMatchesName,saveFlag,f1,d1,f2,d2,TransType)
% SIFT_MOSAIC Demonstrates matching two images using SIFT and RANSAC
%
%   SIFT_MOSAIC demonstrates matching two images based on SIFT
%   features and RANSAC and computing their mosaic.
%
%   SIFT_MOSAIC by itself runs the algorithm on two standard test
%   images. Use SIFT_MOSAIC(IM1,IM2) to compute the mosaic of two
%   custom images IM1 and IM2.

% AUTORIGHTS

% --------------------------------------------------------------------
%                                                         SIFT matches
% --------------------------------------------------------------------

[matches, scores] = vl_ubcmatch(d1,d2) ;
X1 = f1(1:2,matches(1,:)) ; 
X2 = f2(1:2,matches(2,:)) ; 
%check for duplicate matches
%check X1 first
[uX1, IA, IC1] = unique(round(X1'),'rows');
X1 = X1(1:2,IA');
X2 = X2(1:2,IA');
matches = matches(:,IA');

%check X2 next
[uX2, IA, IC1] = unique(round(X2'),'rows');
X1 = X1(1:2,IA');
X2 = X2(1:2,IA');
matches = matches(:,IA');

X1(3,:) = 1 ;
X2(3,:) = 1 ;

%numMatches = size(matches,2) ;
numMatches = size(X1,2);

% --------------------------------------------------------------------
%                                         RANSAC with homography model
% --------------------------------------------------------------------

clear H score ok ;
bestScore = 0;
bestH = zeros(3,3);
bestOK = [];
for t = 1:5000
  % estimate model
  subset = vl_colsubset(1:numMatches, 3) ;
    
  H = eye(3,3);
  if (TransType == 0)
      %Score Using Only Translation
      C1 = mean(X1(1:2,subset),2);
      C2 = mean(X2(1:2,subset),2);
      trans = C2 - C1;
      H(1,3) = trans(1);
      H(2,3) = trans(2);

  elseif (TransType == 1)
      %Score Using Only rotation + translation
      [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', false, 'reflection', false); 
      H = [trans.T' mean(trans.c,1)'; 0 0 1];  

  elseif (TransType == 2)
      %Score Using rotation + translation + scalling
      [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', true, 'reflection', false); 
      H = [trans.b*trans.T' mean(trans.c)'; 0 0 1];  
  
  elseif(TransType == 3)
      %Score using homography
       A = [] ;
       for i = subset
          A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
       end
       [U,S,V] = svd(A) ;
       H = reshape(V(:,9),3,3) ;
  elseif(TransType == 4)
       
      [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', false, 'reflection', false); 
      H = [trans.T' mean(trans.c)'; 0 0 1];
      a=X2(1,subset)*pinv(Z(:,1)');
      b=X2(2,subset)*pinv(Z(:,2)');
      H = [a 0 0; 0 b 0; 0 0 1] * H;
  end
   
   
  X2_ = H * X1 ;
  du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
  dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
  ok = (du.*du + dv.*dv) < 6*6 ;
  score = sum(ok) ;

  if(TransType == 4)
      a=X2(1,ok)*pinv(X2_(1,ok));
      b=X2(2,ok)*pinv(X2_(2,ok));
      H = [a 0 0; 0 b 0; 0 0 1] * H;

      X2_ = H * X1 ;
      du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
      dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
      ok = (du.*du + dv.*dv) < 6*6 ;
      score = sum(ok) ;
  end
  
  if(score > bestScore)
    bestScore = score;
    bestH = H;
    bestOK = ok;  
  end

end

numOkMatches = sum(bestOK);

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

% dh1 = max(size(im2,1)-size(im1,1),0) ;
% dh2 = max(size(im1,1)-size(im2,1),0) ;
% 
% figID=figure(1) ; clf ;
% subplot(2,1,1) ;
% imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]) ;
% o = size(im1,2) ;
% line([f1(1,matches(1,:));f2(1,matches(2,:))+o], ...
%      [f1(2,matches(1,:));f2(2,matches(2,:))]) ;
% title(sprintf('%d tentative matches', numMatches)) ;
% axis image off ;
% 
% subplot(2,1,2) ;
% imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]) ;
% o = size(im1,2) ;
% line([f1(1,matches(1,bestOK));f2(1,matches(2,bestOK))+o], ...
%      [f1(2,matches(1,bestOK));f2(2,matches(2,bestOK))]) ;
% title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
%               sum(bestOK), ...
%               100*sum(bestOK)/numMatches, ...
%               numMatches)) ;
% axis image off ;
% drawnow ;
% colormap('gray')
% if(saveFlag)
%     saveas(figID,saveMatchesName)
% end

% --------------------------------------------------------------------
%                                                               Mosaic
% --------------------------------------------------------------------
% 
% box2 = [1  size(im2,2) size(im2,2)  1 ;
%         1  1           size(im2,1)  size(im2,1) ;
%         1  1           1            1 ] ;
% box2_ = inv(H) * box2 ;
% box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
% box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
% ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
% vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;
% 
% [u,v] = meshgrid(ur,vr) ;
% im1_ = vl_imwbackward(im2double(im1),u,v) ;
% 
% z_ = H(3,1) * u + H(3,2) * v + H(3,3);
% u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
% v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
% im2_ = vl_imwbackward(im2double(im2),u_,v_) ;
% 
% mass = ~isnan(im1_) + ~isnan(im2_) ;
% im1_(isnan(im1_)) = 0 ;
% im2_(isnan(im2_)) = 0 ;
% mosaic = (im1_ + im2_) ./ mass ;
% 
% figure(2) ; clf ;
% imagesc(mosaic) ; axis image off ;
% title('Mosaic') ;
% 
% if nargout == 0, clear mosaic ; end

end