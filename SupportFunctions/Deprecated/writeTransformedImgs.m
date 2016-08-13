
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
        figure(i)
        imshow(imCombined)
        [pathstr,name,ext] = fileparts(char(imageFilename{n})) ;

        %imwrite(im_,strcat(name,'_aligned','.png'),'alpha',double(~isnan(im_))); 
        n
        pause
    end

    figure(i);clf;
    imshow(imCombined)
    imwrite(imCombined,strcat('ref_',num2str(i),'_combined','.png'),'alpha',double(~isnan(imCombined))); 
end
