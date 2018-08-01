function [ ] = git_version_check( username, reponame, checkversion, varargin )
%[ ] = git_version_check( username, reponame, checkversion, varargin )
%
%   This script pulls down information about the latest release of software
%   from a given username and repo, and automatically checks to see if the
%   version we are running (supplied by the user of this function) is less
%   than the target release version on GitHub.
%   
%   Arguments:
%
%   username- The owner of the repository to check against.
%
%   reponame- The name of the repository to check against.
%
%   checkversion- The version we'll be comparing the GitHub release version
%   to. All versions must be preceeded with a 'v'. 
%       Examples: 'v1.0', 'v1.4.3', 'v0.1'
%
%   Parameters:
%
%   'TargetVersion' - The version to check against on GitHub. If not
%   specified, the script will pull down the latest version.
%
%   'WarningOnly' - If true, then no message box will appear.
%
%   Versions are compared left to right, and only compared for as long as
%   the TargetVersion string. So if the checkversion is 'v1.5.6', and the
%   TargetVersion is 'v1.5', the software will assume the checkversion is
%   higher. Conversely, if the checkversion is 'v1.5' and the TargetVersion
%   is 'v1.5.6', then the checkversion will be considered obsolete.
%
%
% Robert F Cooper 07-27-2018

defaultver = 'latest';
defaultmsg = false;

p = inputParser;
p.CaseSensitive = false;
addRequired(p,'username')
addRequired(p,'reponame')
addRequired(p,'checkversion',@ischar)
addParameter(p,'TargetVersion',defaultver,@ischar)
addParameter(p,'WarningOnly',defaultmsg,@islogical)


parse(p,username,reponame,checkversion,varargin{:})
res = p.Results;

repoapisite = ['https://api.github.com/repos/' res.username '/' res.reponame '/releases/' res.TargetVersion];

try
    readfromgit = webread( repoapisite );

    checkvers = strsplit(checkversion(2:end), '.');
    targetvers = strsplit(readfromgit.tag_name(2:end), '.');
    
    
    for i=1:length(targetvers)
        if i>length(checkvers)
            checkver = 0;
        else
            checkver = str2double(checkvers{i});
        end

        targetver = str2double(targetvers{i});
        
        if targetver>checkver
            msgstr = [{['New version available (' readfromgit.tag_name ')! You are running ' checkversion '.']}...
                      {['Go to: https://www.github.com/repos/' res.username '/' res.reponame '/releases/' res.TargetVersion...
                      ' for the latest version.']} ];
            if ~res.WarningOnly
                msgbox(msgstr,'New version available!')
            end
            warning(cell2mat(msgstr))
            
            break;
        elseif targetver<checkver
            warning('Running a newer version of this code than the target version.');
            break;
        end 
    end
    
    
    
catch me
    if strcmp(me.identifier, 'MATLAB:webservices:HTTP404StatusCodeError')
        warning(['Failed to read: ' repoapisite '. Unable to connect to GitHub, OR working from a copy of the repo that does not yet have a release. YMMV!'])
    else
        rethrow(me);
    end
end

end

