# AOAutomontaging
Automontaging for AO images

The enclosed software was developed for montaging adaptive optics images. 

The code was written and tested in MATLAB version R2015b.

Run AutoAOMontagingGUI.m in MATLAB to start the GUI.

--------GUI Manual--------

--Button Descriptions--

Images Folder - Selects an input images folder containing all images that will be montaged. The GUI will automatically look for all confocal, split detection, and dark field images in the folder using the substrings defined in "Input Settings" (see below). 

Scaling/Position File - Selects a scaling/position Excel file that defines the nominal location of each image in the input image folder.See "..\Demos\Multi-Modal\ExamplePositionFile.xlsx" included with the software for the correct format for the file. (This file is not needed for Canon images.)

Output Folder - Selects an output folder where all the montaged images will be saved. 

Montage - Begins the montaging process. The "Images Folder", "Output Folder" and "Scaling/Position File" (if needed) must be filled before beginning the montage.  

--Options Selection--

Montage Type - Selects if the algorithm is starting a new montage or appending to an already completed montage. If appending, the input folder must contain every image used in the previously completed montage along with the new images you would like to append. The output folder must contain the completed montage you would like to append to.

Transformation Type - Select the type of transformation allowed to be used on the images to construct the montage. (The recommended setting is Rigid.)

--Preference Menu--

Input Settings - This will open up a submenu where you can set the substrings used to search for the confocal, split detection, and dark field images in the input images folder.

Output Settings - For future use.

Device Mode - Sets the device that was used to acquire the AO images. Different devices will have different input requirement. E.g., Canon images contain all necessary information in its file name. Hence, the search substrings and Scaling/Position File are not needed. 

Copyright (C) 2016, Min Chen and Robert F. Cooper
All rights reserved.

