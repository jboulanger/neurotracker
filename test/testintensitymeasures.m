% testintensitymeasure
%
% Simple script to test the intensity measurement.
%
% Jerome Boulanger 2017

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

% modify the stage position calibration [FIXIT]
nt.stageposition = nt.stageposition * 5;
figure(1), clf;
imshowpair(log(nt.pano(1,50,0)+0.01),log(nt.pano(2,50,0)+0.01),'falsecolor','ColorChannels', [1 2 0]);

%% compute track and ratio using default parameters
figure(2), clf;
window = 64;
radius = [20 40];
pvalue = [0.1 0.1];
smoothing = 3;
skipframe = 100;
[X,Y,R,resultdata] = trackneuron( nt, window, radius, pvalue, smoothing, skipframe);
% export a log from resultdata
exportlog(strrep(filepath,'.TIFF','.log'), resultdata);
% save directly result data as a .mat file
save(strrep(filepath,'.TIFF','-track.mat'), 'resultdata');

%% represent the track with a color coded intensity ratio
fiexportlog(strrep(filename,'.tif','.log'), resultdata);gure(3), clf;
colorplot(X,Y,R,1)
grid on
box on
xlabel('X ({\mu}m)')
ylabel('Y ({\mu}m)')
title('Color-coded intensity ratio');
colormap jet
colorbar;
