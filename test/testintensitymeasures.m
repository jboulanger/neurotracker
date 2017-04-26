% testintensitymeasure
%
% Simple script to test the intensity measurement
%
% Jerome Boulanger 2017
addpath('../src')
[filename, folder] = uigetfile('*.tiff;*.tif','Pick a calibration file');
filepath = [folder filename];
nt = neurotrackertiff(filepath);
nt.stageposition = nt.stageposition * 5; % why is that so?
[X,Y,R] = trackneuron( nt );