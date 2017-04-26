[filename, folder] = uigetfile('*.tiff;*.tif','Pick a calibration file');
filepath = [folder filename];
nt = neurotrackertiff(filepath);
nt.stageposition = nt.stageposition * 5; % why is that so?
clf, nt.playmovie(true,10);
[X,Y,R] = trackneuron( nt );