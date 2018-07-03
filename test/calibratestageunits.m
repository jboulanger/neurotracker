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
%% load the file
nt = neurotrackertiff(filepath);
%% reset metadata and calibrate the stage
clf
nt = nt.setpixelsize(13); % reset the pixel size as always wrong
nt.magnification = 12.87; % reset the magnification as the estimated by hand on the image
nt = nt.calibratestageunit(19);
nt.printinfo();
%% adjust the flips and play the movie
clf
nt.flip = [true,false;true,false];
nt.playmovie(true, 19);
%% Show a stiched images of all the frames
clf;
imshowpair(nt.pano(1,10,1),nt.pano(2,10,1),'falsecolor','ColorChannels', [1 2 0]);
title('Check calibration');
xlabel('x ({\mu}m)')
ylabel('y ({\mu}m)')

