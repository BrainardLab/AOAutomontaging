function AOAutomontagingLocalHook
% AOAutomontagingLocalHook
%
% Local hook template, toonfigure things for working on the AOAutomontaging project.
%
% For use with the ToolboxToolbox.  If you copy this FILE into your
% ToolboxToolbox localToolboxHooks directory (by default,
% % /<Home>/MATLAB/localHookFolder) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'AOAutomontaging'}) to set up for
% this project.
%
% This will add subfolders of the project to the path as
% well as define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the following project location and i/o directory 
%locations to match what is true on your computer.

%% Say hello
theProject = 'AOAutomontaging';
fprintf('Running %s local hook\n',theProject);

%% Clear out project prefs to avoid staleness
if (ispref(theProject))
    rmpref(theProject);
end

%% You need to set the directories in this cell

% Path to overall AOAutomontaging project folder on your local machine
projectPath = '/Users/Shared/Matlab/Analysis/AOAutomontaging/';

% Input and output parent directories.  The input directories should exist
% and have data in it.  The sample data used in our paper are available at:
%
% https://figshare.com/s/ecaafb98400b0282eaff
% https://figshare.com/s/2cc83c3cffa0cbbdfd33
% https://figshare.com/s/b548a6989c840fabd832

% If you want to reproduce the figures and results in the paper, download
% each sample data zip file and put the directores with names like
% CS_13212_20160104_OS_Images-DONE into the top level of the input
% directory.
%
% These three variables should point to the three downloaded data
% directories
inputDataDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AllInputData/AOMontagingDataSet';
inputManualDataDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AllInputData/AOManualMontages';
inputOverlapAnalysisDataDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AllInputData/OverlapAnalysisPairs';

% Output of the batch processing will show up here.  We make the directory for
% you if it does not exist.
outputMontageDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AOAutomontagingMontageOutput';

% The analysis of montages will be saved to this output directory
outputAnalysisDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AOAutomontagingAnalysisOutput';

%% Make output directories if necessary
if (~exist(outputMontageDir,'dir'))
    mkdir(outputMontageDir);
end
if (~exist(outputAnalysisDir,'dir'))
    mkdir(outputAnalysisDir);
end

%% Set the preferences for the project
setpref(theProject,'inputDataDir',inputDataDir);
setpref(theProject,'inputManualDataDir',inputManualDataDir);
setpref(theProject,'inputOverlapAnalysisDataDir',inputOverlapAnalysisDataDir);
setpref(theProject,'outputMontageDir',outputMontageDir);
setpref(theProject,'outputAnalysisDir',outputAnalysisDir);

