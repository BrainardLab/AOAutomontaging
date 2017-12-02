function saveTif(imageToSave,outputDir,saveFileName)

%make file
tob = Tiff(fullfile(outputDir,saveFileName),'w');

%# you need to set Photometric before Compression
tob.setTag('Photometric',Tiff.Photometric.MinIsBlack)
tob.setTag('Compression',Tiff.Compression.LZW)

%# tell the program that channel 4 is alpha
tob.setTag('ExtraSamples',Tiff.ExtraSamples.AssociatedAlpha)

%# set additional tags (you may want to use the structure
%# version of this for convenience)
tob.setTag('ImageLength',size(imageToSave,1));
tob.setTag('ImageWidth',size(imageToSave,2));
tob.setTag('BitsPerSample',8);
tob.setTag('RowsPerStrip',16);
tob.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Separate);
tob.setTag('Software','MATLAB')
tob.setTag('SamplesPerPixel',2);

if(size(imageToSave,3) == 1)

imageToSave = repmat(imageToSave(:,:,1),[1,1,2]);
imageToSave(:,:,2) = ~isnan(imageToSave(:,:,1));

end

if(isa(imageToSave,'double') || isa(imageToSave,'single'))
    imageToSave = uint8(round(imageToSave*255));
end

%# write and close the file
tob.write(imageToSave)
tob.close
end