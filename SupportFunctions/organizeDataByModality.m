function [inData, MN] = organizeDataByModality(imageDir, ModalitiesSrchStrings)

Allfiles = dir(fullfile(imageDir,'*.tif'));
Allfiles = {Allfiles.name};

MN = size(ModalitiesSrchStrings,1);
inData = [];

%create multi-modal data structure
%depending on parameter array ModalitiesSrchStrings
counter = 0;
for m = 1:MN
    if (~isempty(ModalitiesSrchStrings{m}))%check it's not empty
        imagesFound = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, ModalitiesSrchStrings{m}))));%search for image of this modality
        if(~isempty(imagesFound))
            inData = [inData; imagesFound];
            counter = counter + 1;
        end
        
    end
end
MN = counter;%only using nonempty identifiers

