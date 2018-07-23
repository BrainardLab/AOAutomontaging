function saveTifDouble(imageToSave,outputDir,saveFileName)

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
tob.setTag('BitsPerSample',32);
tob.setTag('RowsPerStrip',16);
tob.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Separate);
tob.setTag('Software','MATLAB')
tob.setTag('SamplesPerPixel',2);
tob.setTag('SampleFormat',3);

%# write and close the file
tob.write(imageToSave)
tob.close
end