% calibratestagerunits
%
% Check the units of the stage position used in the neurotracker
% The conversion is currently hard coded in neurotrackertiff.
%
%
% Jerome Boulanger 2017
clear all; close all;
addpath('../src')
[filename, folder] = uigetfile('*.tiff;*.tif','Pick a calibration file');
filepath = [folder filename];
%% load the file
nt = neurotrackertiff(filepath);
%% reset metadata and calibrate the stage 
clf
nt = nt.setpixelsize(13); % reset the pixel size as always wrong
nt.magnification = 12.87; % reset the magnification as the estimated by hand on the image
stepframe = 5;
dframe = 2;
method = 'imregcorr'; % used either phasecorrelation or compute_displacement
nt = nt.calibratestageunit(stepframe,method, dframe); % use a frame over 20 to calibrate the stage using motion registration
nt.printinfo();
%% adjust the flips and play the movie
clf
nt.playmovie(true, 20);
%% Show a stiched images of all the frames
clf;
imshowpair(nt.pano(1,50,1),nt.pano(2,50,1),'falsecolor','ColorChannels', [1 2 0]);
title('Check calibration');
xlabel('x ({\mu}m)')
ylabel('y ({\mu}m)')
%% Save the stage calibration to a file in the same folder
nt.savestagecalibration([folder 'stagecalibration.mat']);
%% finish
nt.closetiff();
