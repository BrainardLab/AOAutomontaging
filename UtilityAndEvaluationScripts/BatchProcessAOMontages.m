%Script for batch processing AO montages without the GUI
%Written by Min Chen (minchen1@upenn.edu)


%load High Level Input and Output Folders
imageDirBase = getpref('AOAutomontaging','inputDataDir');
outputDirBase = getpref('AOAutomontaging','outputMontageDir');

imageFileFolders = dir(imageDirBase);
for f = 1:length(imageFileFolders)%Do for each data folder
    imageFileFolder          = imageFileFolders(f).name;
    if strcmp(imageFileFolder(1),'.')
        continue;
    end
    
    for TransType = [0 1] %repeat for Translation Only and Rigid Modes 
        
        switch TransType
            case 0
                outputDirTrans = 'TranslationOnly';
            case 1
                outputDirTrans = 'RotationTranslation';
            case 3
                outputDirTrans = 'Homography';
        end
        
        for m = 1:5 %repeat for different combinations of modalities
            switch m
                case 1
                    ModalitiesSrchStrings = {'confocal'};
                    outputDir = 'Confocal';
                case 2
                    ModalitiesSrchStrings = {'split_det'};
                    outputDir = 'SplitDetection';
                case 3
                    ModalitiesSrchStrings = {'avg'};
                    outputDir = 'DarkField';
                case 4
                    ModalitiesSrchStrings = {'confocal'; 'split_det'};
                    outputDir = 'Confocal_and_SD';
                case 5
                    ModalitiesSrchStrings = {'confocal'; 'split_det'; 'avg'};
                    outputDir = 'AllThree';
            end
            
            imageDir=fullfile(imageDirBase,imageFileFolder,'De-identified Data');
            posFileLoc=fullfile(imageDirBase,imageFileFolder,'Seed_Locs_Convert.xls');
            outputDirFull = fullfile(outputDirBase,imageFileFolder,outputDirTrans,outputDir);
            if (~exist(outputDirFull,'dir'))
                mkdir(outputDirFull);
            end
            AOMosiacAllMultiModal(imageDir,posFileLoc,outputDirFull,'multi_modal',ModalitiesSrchStrings,TransType,0,[])
        end
        
    end
end


