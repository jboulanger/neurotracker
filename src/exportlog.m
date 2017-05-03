function exportlog(filename, data)
% exportlog(filename, data)
%
% Export the tracking data as a log file (Zoltan's log format)
%
% Jerome Boulanger 2017
%

h = fopen(filename,'w');
fprintf(h, 'frame bg1 bg2 x1 y1  area1 mean1   x2 y2  area2 mean2\n');
for i=1:length(data)
    fprintf(h, '%d ', i); %frame
    fprintf(h, '%.0f ', data(i).background1.MeanIntensity); %left background
    fprintf(h, '%.0f ', data(i).background2.MeanIntensity); %right background
    fprintf(h, '%.2f ', data(i).signal1.position(1)); %left x
    fprintf(h, '%.2f ', data(i).signal1.position(2)); %left y
    fprintf(h, ' ');
    fprintf(h, '%.0f ', data(i).signal1.Area); %left area
    fprintf(h, '%.0f ', data(i).signal1.MeanIntensity); %left value
    fprintf(h, ' ');
    fprintf(h, ' ');
    fprintf(h, '%.2f ', data(i).signal2.position(1)); %left x
    fprintf(h, '%.2f ', data(i).signal2.position(2)); %left y
    fprintf(h, ' ');
    fprintf(h, '%.0f ', data(i).signal2.Area); %left area
    fprintf(h, '%.0f ', data(i).signal2.MeanIntensity); %left value    
    fprintf(h, '\n');
end
fclose(h);