% calibratestagerunits
%
% Check the units of the stage position used in the neurotracker
% The conversion is currently hard coded in neurotrackertiff.
%
%
% Jerome Boulanger 2017

addpath('../src')
[filename, folder] = uigetfile('*.tiff;*.tif','Pick a calibration file');
filepath = [folder filename];
nt = neurotrackertiff(filepath);
nt.stitch(50)
title('Check calibration with 5cm ruler');
xlabel('x ({\mu}m)')
ylabel('y ({\mu}m)')
