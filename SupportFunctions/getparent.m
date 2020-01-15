function [ parent kids ] = getparent( path, height, returntype)
% Robert Cooper 08-29-11
%   This function returns the parents of a path, up to a height designated
%   by height. It then returns the parent directories in the parent
%   variable and the children in the kids var.
% 
% Input of the return type is useful when trying to get only the parent
% directory. Use 'short' to get only the name of the parent directory of
% height h, or use 'full' to get the entirety of the parent path.
% If the height is undefined it will be assumed to be 1.

if ~exist('height','var') % If they don't input height
    height=1;
    returntype='full';
elseif ischar(height)
    returntype=height;
	height=1;
end


if strcmp(returntype,'short')
    indices=regexp(path,filesep);
    
    if height==0
        parent=path(indices(end)+1:end );
    elseif (height>0) && (height<length(indices))
        parent=path(indices(end-height)+1:indices(end-(height-1))-1 );
    else
        error('Incorrect height value- number must be positive and less than the length of the path');
    end
    
    kids='';
elseif strcmp(returntype,'full')
    indices=regexp(path,filesep);
    
    if height==0 % Included for consistency... no idea why you'd want this though.
        parent=path;
        kids='';
    elseif (height>0) && (height<length(indices))
        parent=path(1:indices(end-(height-1))-1 );
        kids=path(indices(end-(height-1))+1:end );
    else
        error('Incorrect height value- number must be positive and less than the length of the path');
    end
           
end



end

