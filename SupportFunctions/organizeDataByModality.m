function [inData, MN, errorFlag] = organizeDataByModality(imageDir, ModalitiesSrchStrings)

Allfiles = dir(fullfile(imageDir,'*.tif'));
Allfiles = {Allfiles.name};

MN = size(ModalitiesSrchStrings,1);
inData = [];
errorFlag=0;

%create multi-modal data structure
%depending on parameter array ModalitiesSrchStrings
counter = 0;
for m = 1:MN
    if (~isempty(ModalitiesSrchStrings{m}))%check it's not empty
        imagesFound = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, ModalitiesSrchStrings{m}))));%search for image of this modality
        if(~isempty(imagesFound))
            if size(imagesFound,2) == size(inData,2) || isempty(inData)
                inData = [inData; imagesFound];
                counter = counter + 1;
            else 
                errorFlag = ['Error: Mismatch detected. Dataset must have the same number of modalities.'];
                return
            end
        end
        
    end
end
MN = counter;%only using nonempty identifiers

% Go through each column, stripping out modality information. What
% remains should be the same if we have 1:1 correspondence. If not,
% throw an error.
for n = 1:size(inData,2)
    compstrs = cell(MN,1);
    compstrsInd=zeros(MN,1);
    for m=1:MN
        if ~isempty(ModalitiesSrchStrings{m})
            compstrsInd(m) = strfind(inData{m,n},ModalitiesSrchStrings{m});
            compstrs{m} = inData{m,n}(:,1:(compstrsInd(m)-1));
        end
    end
     if ~all(compstrsInd(1) == compstrsInd)
         errorFlag = ['Error: All image pairs/tuples must have identical prefixs before modality identifier. Check image ' inData{1,n} '.'];
         return
     end
    compstrs=cell2mat(compstrs);
    for m=1:MN
        if ~all(compstrs(1,:) == compstrs(m,:))
         errorFlag = ['Error: All image pairs/tuples must have identical prefixs before modality identifier. Check image ' inData{1,n} '.'];
            return
        end
    end
end



