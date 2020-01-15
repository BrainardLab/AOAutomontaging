function cellsizes = findcellsizes(ind,im)

w = 25;%window size left and right
cellsizes = -ones(size(ind));
for i = 1:size(ind)

    [y, x] = ind2sub(size(im),ind(i));
    
    Yrange = max(y-w,1) :1:min(y+w,size(im,1));
    Xrange = max(x-w,1) :1:min(x+w,size(im,2));
    window = im(Yrange,Xrange);
    
    cy = min(w+1, y);
    cx = min(w+1, x);

    window(cy,cx) = 0;
    
    edgedists = bwdist(window);
    
    cellsizes(i) = edgedists(cy,cx);    
    
end


