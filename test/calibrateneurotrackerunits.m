% calibrateneurotrackerunits
%
% Check the units of the stage position used in the neurotracker
%
% Jerome Boulanger 2017

[filename, folder] = uigetfile('*.tiff;*.tif','Pick a calibration file');
filepath = [folder filename];
nt = neurotrackertiff(filepath);
nt.stitch(50)
title('Check calibration with 5cm ruler');
xlabel('x ({\mu}m)')
ylabel('y ({\mu}m)')
