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

end

