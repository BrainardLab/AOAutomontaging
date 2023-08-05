function outNameList = AOMontageLongitudinal(imageDirBase, posFileLocBase, imageDirFollow, posFileLocFollow,outputDir,deviceMode,ModalitiesSrchStrings,...
    BaseTransType,BaseFeatureType,FollowTransType,FollowFeatureType,ExistingMontageBase)


%Montage Baseline
outputDirBase= fullfile(outputDir,'Baseline');
mkdir(outputDirBase);
outputDirFollow= fullfile(outputDir,'Followup');
mkdir(outputDirFollow);

AppendToExisting = 0;
MontageSave = [];
export_to_pshop = 0;

%calculate BaseLine Montage if Needed
if(~exist('ExistingMontageBase','var') || isempty(ExistingMontageBase))
[temp, imNameAll_A, TransformsA, f_all_A, d_all_A] = AOMosiacAllMultiModal(imageDirBase, posFileLocBase, outputDirBase, deviceMode,...
ModalitiesSrchStrings,BaseTransType,AppendToExisting,MontageSave,export_to_pshop,BaseFeatureType);
ExistingMontageBase=fullfile(outputDirBase,'AOMontageSave.mat');
end

[outNameList, TotalTransform, MatchedTo_Followup_to_Baseline] = AOAlignToExistingMontage(imageDirBase, posFileLocBase, ExistingMontageBase, imageDirFollow,...
    posFileLocFollow, outputDirFollow, deviceMode, ModalitiesSrchStrings,FollowTransType,FollowFeatureType)
end