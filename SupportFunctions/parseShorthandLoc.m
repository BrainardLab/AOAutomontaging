function LocXY = parseShorthandLoc(inputString,eyeSide)
%Parses shorthand strings to determine coordinate locations.
%Inputs:
%inputString -- string to parse
%eyeSide -- string indicating which eye is being observed (OD or OS)
%Output:
%LocXY -- [X, Y] location parsed from string
%Written by Min Chen (minchen1@upenn.edu)

inputString = upper(inputString);%make input all caps just in case
LocXY = nan(2,1);

%Flip x coordinate if eye is 'OD'
EyeFlip = 1;
if(strcmpi(eyeSide,'OD'))
    EyeFlip = -1;
end

%check foveal cases first
switch inputString
    
    case 'TL'
        LocXY = [1, 1]';
    case 'TM'
        LocXY = [0, 1]';
    case 'TR'
        LocXY = [-1, 1]';
    case 'ML'
        LocXY = [1, 0]';
    case 'C'
        LocXY = [0, 0]';
    case 'MR'
        LocXY = [-1, 0]';
    case 'BL'
        LocXY = [1, -1]';
    case 'BM'
        LocXY = [0, -1]';
    case 'BR'
        LocXY = [-1, -1]';
end

%if not foveal case then parse
if(isnan(LocXY(1)) || isnan(LocXY(2)))
    
    lettersLoc = find(isletter(inputString));
    LN = length(lettersLoc);
    N = length(inputString);
    if( LN > 2 || LN < 1 || N < 2 || (LN == 2 && N < 4)) %check that format fits 1 or 2 letters each followed by nubmers
        return;
    else%if formats fits then we initialize with zeros
        LocXY = [0 0]';
    end
    
    if(lettersLoc(1) == 1)%Letter infront format
        for i = 1:LN
            
            if(LN == 2 && i == 1)
                lastNum=lettersLoc(2)-1;
            else
                lastNum = N;
            end
            
            switch upper(inputString(lettersLoc(i)))
                case 'T'
                    LocXY(1)=EyeFlip*str2double(inputString(lettersLoc(i)+1:lastNum));
                case 'N'
                    LocXY(1)=-EyeFlip*str2double(inputString(lettersLoc(i)+1:lastNum));
                case 'S'
                    LocXY(2)=str2double(inputString(lettersLoc(i)+1:lastNum));
                case 'I'
                    LocXY(2)=-str2double(inputString(lettersLoc(i)+1:lastNum));
            end
        end
    elseif(lettersLoc(LN) == N)%Letter Behind Format
        for i = 1:LN
            
            if(LN == 2 && i == 2)
                firstNum=lettersLoc(1)+1;
            else
                firstNum = 1;
            end
            
            switch upper(inputString(lettersLoc(i)))
                case 'T'
                    LocXY(1)=EyeFlip*str2double(inputString(firstNum:(lettersLoc(i)-1)));
                case 'N'
                    LocXY(1)=-EyeFlip*str2double(inputString(firstNum:(lettersLoc(i)-1)));
                case 'S'
                    LocXY(2)=str2double(inputString(firstNum:(lettersLoc(i)-1)));
                case 'I'
                    LocXY(2)=-str2double(inputString(firstNum:(lettersLoc(i)-1)));
            end
        end
    end
    
    
    
end


end