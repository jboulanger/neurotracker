# neurotracker
Measure of neuron activity in C. Elegans worm

- Calibrate the stage using the script 'calibratestageunits.m'. It
will determine the camera orientation and the stage calibration by
analysing the motion in an image sequence where a graticule is scanned
in X and Y.
```
nt = nt.setpixelsize(13);
nt.magnification = 12.87;
stepframe = 5;
dframe = 2;
method = 'imregcorr';
nt = nt.calibratestageunit(stepframe,method, dframe);
```
Check the calibration by displaying a montage of the
sequence using the function 'neurotrackertiff.pano()':
```
imshowpair(nt.pano(1,50,1),nt.pano(2,50,1),'falsecolor','ColorChannels', [1 2 0]);
```
Then save the calibration in a file:
```
nt.savestagecalibration([folder 'stagecalibration.mat']);
```

- Analyse intensities and tracks using the script 'testintensitymeasures.m'.
Load the file and the calibration using:
```
nt = neurotrackertiff(filepath);
nt = nt.setpixelsize(13);
nt = nt.setmagnification(12.8);
nt = nt.loadstagecalibration([folder 'stagecalibration.mat']);
```
Then track the neuron and store the intensties in a circular region around the spots of interest:
```
tracks = trackneuron( nt, window, radius, pvalue, smoothing, skipframe);
```
Display the tracks and ratio manually:
```
A = 0.67;
B = 0.93;
[X,Y] = tracks.position();
[R,T] = tracks.ratio(A,B,0);
subplot(121), colorplot(X,Y,R,1)
subplot(122), plot(T,R)
```

## Detail of classes in src/
### neurotrackertiff

Load images from a TIFF produced by neurotracker tiff

### trackset

Store tracks of several spots across the time sequence. Also stores

- To retreive the position over time of all spots:
```
[x,y] = tracks.position();
plot(x(1,:),y(1,:));
```

- To get the ratio with the correction formula (R-A)*B
```
A = 0.67;
B = 0.93;
[r,t] = tracks.ratio(A,B);
plot(r(1,:),t(1,:));
```

## TIFF metadata
- Tag 65000    TAG_TRACKING_POSITION_X,  TIFF_SHORT,  "Tracking Result X component [pixel]" ,
- Tag 65001    TAG_TRACKING_POSITION_Y, TIFF_SHORT,  "Tracking Result Y component [pixel]",
- Tag 65002    TAG_TARGET_POSITION_X,  TIFF_SHORT,  "Target Position X component [pixel]",
- Tag 65003    TAG_TARGET_POSITION_Y, TIFF_SHORT,  "Target Position Y component [pixel]",
- Tag 65004    TAG_TRAY_POSITION_X,  TIFF_LONG,  "Stage Position X component",
- Tag 65005    TAG_TRAY_POSITION_Y, TIFF_LONG,  "Stage Position Y component",
- Tag 65006    TAG_PIXEL_WIDTH,  TIFF_FLOAT,  "Pixel width [micrometer]",
- Tag 65007    TAG_PIXEL_HEIGHT, TIFF_FLOAT,  "Pixel height [micrometer]",
- Tag 65008    TAG_MAGNIFICATION,  TIFF_SHORT,  "Magnification",
- Tag 65009    TAG_TIMESTAMP, TIFF_LONG,  "Timestamp starting from 0 [ms]",
- Tag 65010    TAG_EXPOSURE_TIME,  TIFF_LONG,  "Exposure time [ms]"
