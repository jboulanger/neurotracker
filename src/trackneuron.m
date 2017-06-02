function results = trackneuron( imgsrc, window, radius, pvalue, smoothing, skipframe)
% tracks = TRACKNEURON (imgsrc, window, radius, pvalue, smoothing, skipframe )
% tracks = TRACKNEURON (imgsrc, window, radius, pvalue, smoothing )
% tracks = TRACKNEURON (imgsrc, window, radius, pvalue )
% tracks = TRACKNEURON (imgsrc, window, radius )
% tracks = TRACKNEURON (imgsrc, window )
% tracks = TRACKNEURON (imgsrc )
%
%  Track a neuron from a sequence of images and measure its activity
%
%  The input parameters are:
%  * 'imgsrc' is an object with imread tostagecoords and elapsedtime methods.
%  Meaning that the functions
%   imgsrc.imread(frame,channel)
%   imgsrc.tostagecoords(x,y,frame,channel)
%   imgsrc.elapsedtime(frame, channel)
%  are defined. 
%  You can use imagsrc = neurotrackertiff('path to file') for example.
%  * 'window' is a croping window to accelerate processing
%  * 'radius' is a vector [r1,r2] defining the disks and ring for
%  signal/background intensity analysis
%  * 'pvalue' is a vector [p1,p2] used to detect signal/background intensity
%  * 'smoothing' is the standard deviation of the gaussian filter used to
%  pre-filter the image for detection/segmentation.
%  * 'skipframe' control how often the frames are displayed. skipframe=1
%  shows all frame, skipframe=10, 1 every 10. skipframe<0, disable
%  display.
%
%  The 3 but last parameters are used to feed the function 'measureintensity'
%
%  The output is coordinates lists for each frame X,Y and the intensity ratio
%
%   Jerome Boulanger 2016-2017

if nargin < 2
    window = 64;
end
if nargin < 3
    radius = [20, 40];
end
if nargin < 4
    pvalue = [0.1, 0.1];
end
if nargin < 5
    smoothing = 3;
end
if nargin < 6
    skipframe = 1;
end

results = trackset(1,imgsrc.length(),2);

tic;
for t=1:imgsrc.length(); 
    im1 = double(imgsrc.imread(t,1));
    w = size(im1,1);
    if window > 0
        cim1 = im1(w/2-window:w/2+window-1,w/2-window:w/2+window-1);    
    else
        cim1 = im1;
    end
    try 
        [x, y] = detectspots(cim1, smoothing, 1e-3, true);    
    catch 
        x = X(t-1);
        y = Y(t-1);
    end
    [v1,b1] = measureintensity(cim1,x,y,radius,pvalue,smoothing);
    
    im2 = double(imgsrc.imread(t,2));
    if window > 0
        cim2 = im2(w/2-window:w/2+window-1,w/2-window:w/2+window-1); 
    else
        cim2 = im2;
    end
    [v2,b2] = measureintensity(cim2,x,y,radius,pvalue,smoothing);
        
    x = x + w/2 - window;
    y = y + w/2 - window;
    [X(t),Y(t)] = imgsrc.tostagecoords(x,y,t,1);      
    
    R(t) = (v1.MeanIntensity - b1.MeanIntensity) / (v2.MeanIntensity - b2.MeanIntensity);
           
    results = results.set(1, t, 1, makerecord(v1,b1,imgsrc,t,1));
    results = results.set(1, t, 2, makerecord(v2,b2,imgsrc,t,1));
    
    if mod(t,skipframe) == 0 && skipframe > 0
        showchannels(cim1,v1,b1,1,'channel 1');
        showchannels(cim2,v2,b2,2,'channel 2');
        showgraphs(imgsrc, im1, im2, t, X, Y,x,y,radius, R);
    end
end
toc

end
function showchannels(im,v,b,pos,str) 
    subplot(2,2,pos)
    maskv = zeros(size(im));
    maskv(v.PixelIdxList) = im(v.PixelIdxList); 
    maskb = zeros(size(im));
    maskb(b.PixelIdxList) = im(b.PixelIdxList);
    imshowpair(maskv,maskb,'falsecolor','ColorChannels', [1 2 0]);            
    hold on
    plot(0,0,'r.');
    plot(0,0,'g.');
    hold off
    title(str);
    legend('signal','background')
    
end

function showgraphs(imgsrc, im1, im2, t, X, Y,x,y,radius,R)
    subplot(223)
    ref1 = imgsrc.imref2d(t,1);
    ref2 = imgsrc.imref2d(t,2);
    imshowpair(im1,ref1,im2,ref2,'falsecolor','ColorChannels', [1 2 0]);
    hold on
    plot(imgsrc.stageposition(1:2:end,1), imgsrc.stageposition(1:2:end,2));    
    plot(X(1:t),Y(1:t),'g-');
    ex1 = x+radius(1)*cosd(0:360);
    ex2 = x+radius(2)*cosd(0:360);
    ey1 = y+radius(1)*sind(0:360);
    ey2 = y+radius(2)*sind(0:360);
    [ex1, ey1] = imgsrc.tostagecoords(ex1,ey1,t,1);
    [ex2, ey2] = imgsrc.tostagecoords(ex2,ey2,t,1);
    plot(ex1,ey1,'m');
    plot(ex2,ey2,'m--');    
    hold off    
    axis([X(t)-500,X(t)+500,Y(t)-500,Y(t)+500]);     
    xlabel('X ({\mu}m)')
    ylabel('Y ({\mu}m)')
    grid on
    box on
    title('Centered view in stage referential')
    legend('stage','track','signal','background')
    subplot(224)
    plot(imgsrc.elapsedtime(1:t,1)/1000, R(1:t));
    axis([0, imgsrc.elapsedtime(imgsrc.length(),1)/1000, 0, 10])
    xlabel('Time (s)')
    ylabel('ch1/ch2');
    grid on
    box on
    title('Intensity ratio')
    drawnow  
end

function rec = makerecord(I,B,imgsrc,t,c)
    x = I.WeightedCentroid(2);
    y = I.WeightedCentroid(1);   
    [x,y] = imgsrc.tostagecoords(x,y,t,c);
    rec.signal = I;
    rec.background = B;    
    rec.x = x;
    rec.y = y;
    rec.t = imgsrc.elapsedtime(t,c) / 1000;
end

