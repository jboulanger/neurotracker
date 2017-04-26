% testintensitymeasure
%
% Simple script to test the intensity measurement
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

% create a neurotrackertiff object
nt = neurotrackertiff(filepath);
% modify the stage position calibration [FIXIT]
nt.stageposition = nt.stageposition * 5;

% compute track and ratio using default parameters
figure(1)
[X,Y,R] = trackneuron( nt );

% represent the track with a color coded intensity ratio
figure(2)
colorplot(X,Y,R,1)
grid on
box on
xlabel('X ({\mu}m)')
ylabel('Y ({\mu}m)')
title('Color-coded intensity ratio');
colormap jet
colorbar;
