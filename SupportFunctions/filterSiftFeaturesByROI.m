function [f_out, d_out] = filterSiftFeaturesByROI(im, f, d, ROICropPct)
%Filters SIFT features (f and d) to only those in the center (1-2*ROICropPct) of the image.

[rows, cols, planes] = size(im);
xmin=fix(cols*ROICropPct);
ymin=fix(rows*ROICropPct);
xmax=fix(cols*(1-ROICropPct));
ymax=fix(rows*(1-ROICropPct));
Ind=find(f(1,:) >= xmin & f(1,:) <= xmax & f(2,:) >= ymin & f(2,:) <= ymax);


f_out = f(:,Ind);
d_out = d(:,Ind);
cleanEdgesFlag = 0;
if(cleanEdgesFlag)
%finds biggest connected object, erode 50 pixels to the side
%then removes any features in those edges.
biggest = bwareafilt(im>0,1);
ROI = imerode(biggest,ones(50,50));
ROIatLocation = zeros(1,size(f_out,2));
%for each detected feature, we check if it falls inside ROI
for i = 1:length(ROIatLocation)
    ROIatLocation(i) = ROI(round(f_out(2,i)),round(f_out(1,i)));
end
%only keep those inside    
Ind = find(ROIatLocation);
f_out = f_out(:,Ind);
d_out = d_out(:,Ind);
end

end

