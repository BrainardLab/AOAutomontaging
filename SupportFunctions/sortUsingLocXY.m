function [LocXY, inData, imageFilename, pixelScale,ID] = sortUsingLocXY(LocXY, inData, imageFilename, pixelScale, ID, MN)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[~, I] = sortrows(abs(LocXY)',[1,2]);
LocXY=LocXY(:,I);
ID = ID(I);
pixelScale= pixelScale(I);

for m = 1:MN
    imageFilename(m,:)=imageFilename(m,I);
    inData(m,:)=inData(m,I);
end

end

