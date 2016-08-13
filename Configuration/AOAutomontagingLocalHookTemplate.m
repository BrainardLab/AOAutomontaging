function AOAutomontaging
% AOAutomontaging
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'AOMontaging'}) to set up for
% this project.  You then edit your local copy to match your 
%
% The thing that this does is add subfolders of the project to the path as
% well as define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.
%
% Configure things for working on the IBIOColorDetect project.

%% Say hello
fprintf('Running AOAutomontaging local hook\n');

%% You need to set the directories in this cell

% Path to AOAutomontaging project on your local machine
projectPath = '/Users/Shared/Matlab/Analysis/AOAutomontaging/';

% Input and output parent directories.  The input directories should exist
% and have data in it.  Sample data used in our paper is available at
%   https://figshare.com/s/ecaafb98400b0282eaff
%
% If you want to reproduce the figures and results in the paper, downlaod
% the sample data zip file and put the directores with names like
% CS_13212_20160104_OS_Images-DONE into the top level of the input
% directory.
%
% inputDataDir = 'C:\Users\dontm\Documents\Research\AdaptiveOpticsMosaic\PaperValidationExperiments\BOE_2016\Data';
inputDataDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AOAutomontagingInput';

% Output of batch processing will show up here.  We make the directory for
% you if it does not exist.
%
% outputMontageDir = 'C:\Users\dontm\Documents\Research\AdaptiveOpticsMosaic\PaperValidationExperiments\BOE_2016\AutoMontageResults';
outputMontageDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AOAutomontagingMontageOutput';

% Analysis of montages output directory
%
% outputAnalysisDir = 'C:\Users\dontm\Documents\Research\AdaptiveOpticsMosaic\PaperValidationExperiments\BOE_2016\Analysis;
outputAnalysisDir = '/Volumes/Users1/DropboxLab/AOSLOImageProcessing/AOAutomontagingAnalysisOutput';

%% Put project toolbox onto path
tbDeployToolboxes('config',tbToolboxRecord( ...
    'name', 'AOMontagingSupportFunctions', ...
    'type', 'local', ...
    'url', fullfile(projectPath,'SupportFunctions')) ...
    );

%% Make output directories if necessary
if (~exist(outputMontageDir,'dir'))
    mkdir(outputMontageDir);
end
if (~exist(outputAnalysisDir,'dir'))
    mkdir(outputAnalysisDir);
end

%% Set the preferences for the project
setpref('AOAutomontaging','inputDataDir',inputDataDir);
setpref('AOAutomontaging','outputMontageDir',outputMontageDir);
setpref('AOAutomontaging','outputAnalysisDir',outputAnalysisDir);

