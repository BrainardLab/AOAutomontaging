function saveConeLoc(SaveName,CNNPos,imageSize)
%This function allows cone locations to be saved in a parloop

 save(SaveName,'CNNPos','imageSize');