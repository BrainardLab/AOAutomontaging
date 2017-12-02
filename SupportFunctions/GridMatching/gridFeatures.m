function [f,d] = gridFeatures(im)


[im1_mask, X1, Y1] =  findlocalmin(im);


I1 = sub2ind(size(im), Y1, X1);
im1_int = im(I1);

I1_filtered =  I1(im1_int > prctile(im1_int,50));%> mean(im1_int)-std(double(im1_int)));
im1_filtered = zeros(size(im));
im1_filtered(I1_filtered) = 1;

cellsizes1 = findcellsizes(I1_filtered,im1_filtered);
I1_filteredsize =  I1_filtered(cellsizes1 > 6);%prctile(cellsizes1,5));%> mean(im1_int)-std(double(im1_int)));
im1_filteredsize = zeros(size(im));
im1_filteredsize(I1_filteredsize) = 1;

[f, d]=findgridFeature(I1_filteredsize,im1_filteredsize,im,im1_filtered,im1_mask);
