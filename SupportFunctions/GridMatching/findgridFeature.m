function [f d d_int thetas]=findgridFeature(ind, im,im_int,im_filt1, im_filt2)

w=50;%window width on one side
g=5;%grid width in pixels
R=3;%number of Nearest points for rotation variants, including itself(zero rotation)

N=size(ind,1)*(R);% number of features
ND = ceil((2*w+1)/g);%Grid size total;

thetas=nan(1,N); %rotation of feature
f=zeros(2,N); %location of feature
d=zeros(ND^2,N); %feature discriptor

gridarrayTemplate = zeros(ND,ND);%grid size
d_int = zeros((2*w+1)^2,N,4);%debug variable

for i = 1:size(ind)
    
    [y, x] = ind2sub(size(im),ind(i));
    if((y-w >= 1) && (y + w <= size(im,1)) && (x-w >= 1) && (x + w <= size(im,2)))%only care about non-edge cases
        Yrange = y-w :1: y+w;
        Xrange = x-w :1: x+w;
        windowArray = im(Yrange,Xrange);
        windowArray_int = im_int(Yrange,Xrange);
        windowArray_filt1 = im_filt1(Yrange,Xrange);
        windowArray_filt2 = im_filt2(Yrange,Xrange);
        
        %find points in window
        [wy wx] = find(windowArray);
        
        cx = (w+1)*ones(size(wx));%center point
        cy = (w+1)*ones(size(wy));%center point
        
        %vectors from center to each point in window
        x1 = cx - wx;
        y1 = cy - wy;
        
        dist = abs(x1) + abs(y1);%dist to each point
        [sortDist I] = sort(dist);%sort dist by ascending order
        
        %use closest R points to create rotated features
        %***Remember that the center point is included in here
        %so there is at least a zero distance
        
        %rotate to unit vector (makes feature rotation invariant)
        %unit vector we align to
        x2 = 1;
        y2 = 0;
        
        %defaults to zero(no rotation)
        %index R+1 is unrotated feature, so is always zero
        
        for r = 1:min(R,size(wx)) %not always R point, hence the 'min'
            i_r = ((i-1)*R)+r; % this is the array index for the current feature
            
            %use r-th closest vector
                x1_r=x1(I(r));
                y1_r=y1(I(r));
            
            if (x1_r == 0 && y1_r == 0) %we don't rotate when the closest pointing is itself 
                theta = 0;
            else %find angle between vectors
                theta = atan2(x1_r*y2-y1_r*x2,x1_r*x2+y1_r*y2);
            end
            thetas(i_r) = theta; %save theta for output here
            
            %create rotation and translation matrix matrix
            rot = eye(3,3);
            transF = eye(3,3);
            transR = eye(3,3);
            
            rot(1,1)=cos(theta);
            rot(1,2)=-sin(theta);
            rot(2,1)=sin(theta);
            rot(2,2)=cos(theta);
            transF(1,3)=-w-1;
            transF(2,3)=-w-1;
            transR(1,3)=w+1;
            transR(2,3)=w+1;
            
            %H is the final matrix, it moves the center to origin, rotates, then returns to the center
            H = transR*rot*transF;
            %apply transform matrix
            W=[wx';wy';ones(size(wx'))];
            W_trans=H*W;
            
            %convert to grid representation.
            gy = ceil(W_trans(2,:)/g);
            gx = ceil(W_trans(1,:)/g);
            
            gridarray = gridarrayTemplate;
            
            for gi = 1:size(gx,2)
                if((gx(gi) >=1) && (gx(gi) <= size(gridarray,2)) && (gy(gi) >=1) && (gy(gi) <= size(gridarray,1)))
                    gridarray(gy(gi),gx(gi)) = 1;%gridarray(gy(gi),gx(gi))+1;
                end
            end
            % gridarray =  imgaussfilt(gridarray,.75);
            
            %save features and location here
            f(1:2,i_r) = [x ; y];
            d(:,i_r) = gridarray(:);
            d = logical(d);
            d_int(:,i_r,1) = windowArray(:);
            d_int(:,i_r,2) = windowArray_int(:);
            d_int(:,i_r,3) = windowArray_filt1(:);
            d_int(:,i_r,4) = windowArray_filt2(:);
        end
    end
end

%now we clean up features where theta is NaN (i.e. points in window < R)

f=f(:,~isnan(thetas)); %location of feature
d=d(:,~isnan(thetas)); %feature discriptor
d_int = d_int(:,~isnan(thetas),:);
thetas=thetas(~isnan(thetas)); %rotation of feature


