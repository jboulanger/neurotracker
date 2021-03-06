% testintensitymeasure
%
% Simple script to test the intensity measurement.
%
% Jerome Boulanger 2017
clear all; close all;
addpath('../src')

% Select input. Select multiple files is the images are stored in more than
% one file.
[filenames, folder] = uigetfile('*.tiff;*.tif','Select a file','MultiSelect', 'on');

% create the list of filepath if needed
if iscell(filenames)
    filepath = cell(1,numel(filenames));
    for i = 1:numel(filenames)
        filepath{i} = [folder filenames{i}];
    end
else
    filepath = [folder filenames];
end

%% create a neurotrackertiff object
nt = neurotrackertiff(filepath);
nt = nt.setpixelsize(13); % Hamamatsu flash 4 [ 6.5um ] with binning 2x2
nt = nt.setmagnification(12.8); % calibration ??
% load calibration and flip parameters
nt = nt.loadstagecalibration([folder 'stagecalibration.mat']);
nt.printinfo();
% check the stage calibration by creating a panoramic view (stiched)
disp('Create a stiched image');
figure(1), clf;
p1 = nt.pano(1,5,0);
p2 = nt.pano(2,5,0);
imshowpair(log(p1-min(p1(:))+.1),log(p2-min(p2(:))+.1),'falsecolor','ColorChannels', [1 2 0]);

title('Stitched image [please check if the stage registration is correct]')

%% compute track and ratio using default parameters
disp('Extract a track');
scrsz = get(groot,'ScreenSize');
f = figure(2);
set(f,'Position',scrsz), clf;
window = 64;
radius = [20 40];
pvalue = [0.05 0.1];
smoothing = 3;
skipframe = 100;
tracks = trackneuron( nt, window, radius, pvalue, smoothing, skipframe);
disp('Export results')
tracks.save(strrep(filepath,'.TIFF','-track.mat'));
tracks.exportlog(strrep(filepath,'.TIFF','-track.log'));

%% represent the track with a color coded intensity ratio
disp('Display the track');
figure
A = 0.67;
B = 0.93;
ratio_smoothing = 3;
[X,Y] = tracks.position();
[R,T] = tracks.ratio(A,B, ratio_smoothing);
subplot(221)
colorplot(X,Y,T,1)
grid on
box on
xlabel('X ({\mu}m)')
ylabel('Y ({\mu}m)')
title('color coded time');
colormap jet
colorbar;

subplot(222)
colorplot(X,Y,R,1)
grid on
box on
xlabel('X ({\mu}m)')
ylabel('Y ({\mu}m)')
title('color coded ratio');
colormap jet
colorbar;

subplot(212)
plot(T,R);
axis tight
grid on;
title('Ratio over time')
xlabel('Time (s)')
ylabel('Ratio')
%% Interactive display of tracks
inspecttracks(tracks);
disp('Done');
%% Interactive display of tracks
inspectbehavior(tracks);
disp('Done');

%% finish
nt.closetiff();

