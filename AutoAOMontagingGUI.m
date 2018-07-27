function varargout = AutoAOMontagingGUI(varargin)
% AUTOAOMONTAGINGGUI MATLAB code for AutoAOMontagingGUI.fig
%      AUTOAOMONTAGINGGUI, by itself, creates a new AUTOAOMONTAGINGGUI or raises the existing
%      singleton*.
%
%      H = AUTOAOMONTAGINGGUI returns the handle to a new AUTOAOMONTAGINGGUI or the handle to
%      the existing singleton*.
%
%      AUTOAOMONTAGINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOAOMONTAGINGGUI.M with the given input arguments.
%
%      AUTOAOMONTAGINGGUI('Property','Value',...) creates a new AUTOAOMONTAGINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoAOMontagingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoAOMontagingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoAOMontagingGUI

% Last Modified by GUIDE v2.5 31-Aug-2017 14:05:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AutoAOMontagingGUI_OpeningFcn, ...
    'gui_OutputFcn',  @AutoAOMontagingGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end



if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AutoAOMontagingGUI is made visible.
function AutoAOMontagingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoAOMontagingGUI (see VARARGIN)

% Choose default command line output for AutoAOMontagingGUI
handles.output = hObject;

%allocate variables and defaults
handles.outputFolder_name=[];
handles.combinedFile_names=[];
handles.postionFile_name=[];
handles.imgfolder_name=[];
handles.imageFile_names=[];
handles.modalitiesInfo = {'Confocal' 'confocal';
    'Split Detection' 'split';
    'Dark Field' 'avg';
    'Modality 4' '';
    'Modality 5' '';};
handles.inputExt = 1;
handles.device_mode = 'multi_modal';
%default to .tif

%add path and setup vl_feat
currentFile = mfilename('fullpath');
[currentFileLoc,name,ext] = fileparts(currentFile); 
genpath(fullfile(currentFileLoc,'SupportFunctions'));
addpath(genpath(fullfile(currentFileLoc,'SupportFunctions')));
vl_setup;

% If we've set up photoshop, then enable and check the export to photoshop
% buttons.
if exist('psnewdoc')
    set(handles.pshop_cbox,'Enable','on');
    set(handles.pshop_cbox,'Value',1.0);    
end

% Check this version of AO Montaging against git.
fid = fopen(fullfile(getparent(which(mfilename)),'.VERSION'),'r');
if fid ~= -1
    thisver = fscanf(fid,'%s');
    fclose(fid);
    
    git_version_check( 'BrainardLab','AOAutomontaging', thisver )
else
    warning('Failed to detect .VERSION file, unable to determine if running the newest version.')
end

%set default options
set(handles.uibuttongroup1,'selectedobject',handles.radiobutton3);
set(handles.radiobutton4,'enable','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoAOMontagingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AutoAOMontagingGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles)
% hObject    handle to filemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in imageList.
function imageList_Callback(hObject, eventdata, handles)
% hObject    handle to imageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageList

if(~isempty(handles.imgfolder_name) && ~isempty(handles.imageFile_names))
    index_selected = get(handles.imageList,'Value');
    axes(handles.canvas);
    img = imread(fullfile(handles.imgfolder_name,handles.imageFile_names{index_selected}));
    imagesc(img); colormap gray; axis equal; axis off;
end

% --- Executes during object creation, after setting all properties.
function imageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in imagefolder.
function imagefolder_Callback(hObject, eventdata, handles)
% hObject    handle to imagefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selectedDir=uigetdir(handles.imgfolder_name);
if(selectedDir == 0)
    return
end
handles.imgfolder_name = selectedDir;
set(handles.selectFolderText, 'String', handles.imgfolder_name) ;
set(handles.selectFolderText, 'TooltipString', handles.imgfolder_name) ;

Allfiles = dir(fullfile(handles.imgfolder_name,'*.tif'));
Allfiles = {Allfiles.name};

handles.imageFile_names =[];
%Use current identifiers to locate all images
dataSummary = cell(1,1);
dataSummary{1} = 'Input Data Summary:';

if strcmp(handles.device_mode, 'multi_modal')
    for m = 1:size(handles.modalitiesInfo,1)
        if (~isempty(handles.modalitiesInfo{m,2}))%check it's not empty
            found = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, handles.modalitiesInfo{m,2}))));
            handles.imageFile_names = [handles.imageFile_names, found];%search
            dataSummary{end+1} =[num2str(size(found,2)),' ',char(handles.modalitiesInfo{m,1}),' image(s) found.'];
        end
    end
    
    dataSummary{end+1} = '';
    dataSummary{end+1} = 'Note: You can modify how to search for each modality under Preferences->Input Settings.';

elseif strcmp(handles.device_mode, 'canon')
    found = sort(Allfiles(cellfun(@(s) strcmp(s(1:4),'206-'), Allfiles )));
    handles.imageFile_names = [handles.imageFile_names, found];%search
    dataSummary{end+1} =[num2str(size(found,2)),' image(s) found.'];
end

msgbox(dataSummary,'Input Complete');


set(handles.imageList,'String',handles.imageFile_names,...
    'Value',1)
guidata(hObject, handles);


% --- Executes on button press in selectPosFile.
function selectPosFile_Callback(hObject, eventdata, handles)
% hObject    handle to selectPosFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
defaultfolder = handles.imgfolder_name;
[FileName,PathName] = uigetfile(fullfile(defaultfolder,'*.xlsx'));
handles.postionFile_name = fullfile(PathName,FileName);
set(handles.posFileText, 'String', handles.postionFile_name);
set(handles.posFileText, 'TooltipString', handles.postionFile_name);
guidata(hObject, handles);

% --- Executes on selection change in montageList.
function montageList_Callback(hObject, eventdata, handles)
% hObject    handle to montageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns montageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from montageList

if(~isempty(handles.combinedFile_names) && ~isempty(handles.outputFolder_name))
    index_selected = get(handles.montageList,'Value');
    axes(handles.canvas);
    img = imread(fullfile(handles.outputFolder_name,handles.combinedFile_names{index_selected}));
    imagesc(img(:,:,1)); colormap gray; axis equal; axis off;
end


% --- Executes during object creation, after setting all properties.
function montageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to montageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in montageAll.
function montageAll_Callback(hObject, eventdata, handles)
% hObject    handle to montageAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get Transformation Type
TransType=0;
switch get(get(handles.uibuttongroup1,'SelectedObject'),'Tag')
    case 'radiobutton2',  TransType=0;
    case 'radiobutton3',  TransType=1;
    case 'radiobutton4',  TransType=3;
end


%New Montage Or Append to Existing?
AppendToExisting=0;
MontageSave=[];
switch get(get(handles.uibuttongroup2,'SelectedObject'),'Tag')
    case 'radiobutton5'
        % If new check for existing save file exists
        if(exist(fullfile(handles.outputFolder_name, 'AOMontageSave.mat'), 'file'))
            choice = questdlg('AOMontageSave.mat file found in output folder. Overwrite existing montage?','Warrning');
            if(~isequal(choice,'Yes'))
                return;
            end
        end
        AppendToExisting=0;
    case 'radiobutton6',
        %If appending check for existing save file
        if(~exist(fullfile(handles.outputFolder_name, 'AOMontageSave.mat'), 'file'))
            
            
            
            
            choice = questdlg('Appending Failed: AOMontageSave.mat file not found in output folder. Continue as New?','Error');
            if(isequal(choice,'Yes'))
                AppendToExisting = 0;
            else
                return;
            end
        else
            MontageSave = fullfile(handles.outputFolder_name, 'AOMontageSave.mat');
            AppendToExisting=1;
        end
end

%read filename substrings to search for different modalities


tic
handles.combinedFile_names = AOMosiacAllMultiModal(handles.imgfolder_name, handles.postionFile_name, ...
                                                   handles.outputFolder_name, handles.device_mode, ...
                                                   handles.modalitiesInfo(:,2), TransType,AppendToExisting, ...
                                                   MontageSave, get(handles.pshop_cbox,'Value') );
toc

if(~isempty(handles.combinedFile_names))
set(handles.montageList,'String',handles.combinedFile_names,...
    'Value',1)
img=imread(fullfile(handles.outputFolder_name,handles.combinedFile_names{1}));
axes(handles.canvas);
imagesc(img(:,:,1)); colormap gray; axis equal; axis off;
guidata(hObject, handles);
else
errordlg('Montage Interrupted!','Error');    
end




% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;


% --- Executes on button press in outputFolder.
function outputFolder_Callback(hObject, eventdata, handles)
% hObject    handle to outputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selectedDir=uigetdir;
if(selectedDir == 0)
    return
end
handles.outputFolder_name = selectedDir;
set(handles.outputFolderText, 'String', handles.outputFolder_name) ;
set(handles.outputFolderText, 'TooltipString', handles.outputFolder_name) ;
guidata(hObject, handles);


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function prefmenu_Callback(hObject, eventdata, handles)
% hObject    handle to prefmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function inputsettings_Callback(hObject, eventdata, handles)
% hObject    handle to inputsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Filename substrings for identifying each modality:','Image extension:'};
dlg_title = 'Input';
num_lines = [5 50;1 10];
%defaultsettings = {char({'confocal';'split';'avg';''});'.tif'};
defaultsettings = {'confocal\nsplit\navg\n','.tif'};

Title = 'Input Settings';

%%%% SETTING DIALOG OPTIONS
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'off';
Options.ButtonNames = {'Save','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 4; % Horizontal dimension in fields

Prompt = {};
Formats = {};
DefAns = struct([]);

Prompt(1,:) = {'Filename substrings for identifying each modality:','modalitiesInfo',[]};
Formats(1,1).type = 'table';
Formats(1,1).format = {'char', 'char'}; % table (= table in main dialog) / window (= table in separate dialog)
Formats(1,1).items = {'Modality Name' 'Substring'};
Formats(1,1).size = [158.5 106];
Formats(1,1).margin = 2;
Formats(1,1).span = [1 1];  % item is 2 field x 1 fields
Formats(1,1).labelloc = 'topcenter';
Formats(1,1).unitsloc = 'bottomcenter';
DefAns(1).modalitiesInfo = handles.modalitiesInfo;

Prompt(end+1,:) = {'Input images extension:','inputExt',[]};
Formats(2,1).type = 'list';
Formats(2,1).style = 'popupmenu';
Formats(2,1).items = {'.tif','.png'};
DefAns(1).inputExt = handles.inputExt;

[Input,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

if(~Cancelled)
    handles.modalitiesInfo=Input.modalitiesInfo;
    handles.inputExt=Input.inputExt;
    guidata(hObject, handles);
end
% --------------------------------------------------------------------
function outputsettings_Callback(hObject, eventdata, handles)
% hObject    handle to outputsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function multi_modal_device_Callback(hObject, eventdata, handles)
% hObject    handle to multi_modal_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.canon_device, 'Checked', 'off');
set(handles.multi_modal_device, 'Checked', 'on');
set(handles.selectPosFile, 'Enable','on');
set(handles.inputsettings,'Enable','on');

handles.modalitiesInfo = {'Confocal' 'confocal';
                          'Split Detection' 'split';
                          'Dark Field' 'avg';
                          'Modality 4' '';
                          'Modality 5' '';};
handles.device_mode = 'multi_modal';
guidata(hObject, handles);

% --------------------------------------------------------------------
function canon_device_Callback(hObject, eventdata, handles)
% hObject    handle to canon_device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.canon_device, 'Checked', 'on');
set(handles.multi_modal_device, 'Checked', 'off');
set(handles.posFileText,'String','');
set(handles.selectPosFile, 'Enable','off');
set(handles.inputsettings,'Enable','off');

handles.modalitiesInfo = {'Canon confocal','206-'};
handles.device_mode = 'canon';
guidata(hObject, handles);


% --- Executes on button press in pshop_cbox.
function pshop_cbox_Callback(hObject, eventdata, handles)
% hObject    handle to pshop_cbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pshop_cbox
