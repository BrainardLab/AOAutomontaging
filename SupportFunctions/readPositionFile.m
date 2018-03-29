function [imageFilename, eyeSide, pixelScale, LocXY, ID, NC, errorFlag] = readPositionFile(imageDir, inData, posFileLoc,device_mode, matchexp, MN, N)

%defaults
imageFilename = cell(MN,N);%stores all filenames
pixelScale = NaN(1,N); % Stores the relative size of our images
ID = cell(1,N);
LocXY = nan(2,N);%set NaN to start
NC = -1;%number of columns in the file
eyeSide = 'OS';
errorFlag = 0;

if strcmp(device_mode, 'multi_modal')
    
    %load position info from excel spreadsheet
    [temp,temp2,C] = xlsread(posFileLoc);
    C(cellfun(@(x) any(isnan(x)), C)) = {''};
    C = cellfun(@num2str, C,'UniformOutput',false); % First make them all strings.
    NC = size(C,2);
    %First look for key info like which eye
    for i = 1:size(C,1)
        if(strcmpi(C{i,1},'eye'))
            eyeSide = C{i,2};
        end
    end
    
    
    % Then convert back to a number, before adding the trappings of our
    % file ids.
    C(:,1) = cellfun(@(x) ['_' num2str( str2double(x),'%04.0f') '_'], C(:,1),'UniformOutput',false);
    
    %verify that the image id's line up for all modalities
    for n = 1:N
        %build filename structure and check to make sure all modalities are
        %present for all ids
        imageFilename{1,n} = fullfile(imageDir, inData{1,n});
        ImageID_m1 = regexpi(inData{1,n},matchexp,'match');
        for m = 2:MN
            imageFilename{m,n} = fullfile(imageDir, inData{m,n});
            ImageID_mf = regexpi(inData{m,n},matchexp,'match');
            if(~strcmpi(ImageID_m1,ImageID_mf))%check for errors
                errorFlag = ['Error: Mismatch detected. Every image number must have the same number of modalities. Check image ' ImageID_m1];
                return;
            end
        end
        %match with info from excel
        for i = 1:size(C,1)
            if (~isempty(strfind(inData{1,n}, C{i,1})))
                
                if(NC >= 4)
                    scale=str2double(strtrim(C{i,4}));
                    
                    
                    if(~isnan(scale) && scale>0)
                        pixelScale(n) = scale;
                    end
                end
                
                %first try looking at coordinate grid
                if(NC >= 3)
                    ID{n} = C{i,2};
                    Loc = strsplit(C{i,3},',');
                    if(size(Loc,2) == 2)
                        LocXY(1,n) = str2double(strtrim(Loc{1}));
                        LocXY(2,n) = str2double(strtrim(Loc{2}));
                    end
                end
                
                if(NC >= 2)
                    %coordinate grind c
                    if(isnan(LocXY(1,n)) || isnan(LocXY(2,n)))
                        LocXY(:,n) = parseShorthandLoc(C{i,2},eyeSide);
                    end
                end

                break;
            end
        end
        % If we can't find
        if any(isnan(LocXY(:,n))) || isnan(pixelScale(n))
            warning(['Error: Location information missing or invalid for image cluster: ' imageFilename{1,n}]);
            for m = 1:MN
                imageFilename{m,n} = [];
            end
        end
        
    end
    
elseif strcmp(device_mode, 'canon')
    for n = 1:N
        [ eyeSide, LocXY(:,n) ] = parseCanonFName( inData{1,n} );
        imageFilename{1,n} = fullfile(imageDir, inData{1,n});
    end
end

pixelScale = max(pixelScale)./pixelScale;