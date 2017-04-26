function [X,Y,R] = trackneuron( imgsrc )
% tracks = TRACKNEURON (imgsrc)
%   Track a neuron from a sequence of images and measure its activity
%
%   Output: coordinates X,Y  and ratio
%   Jerome Boulanger 2016-2017

X = zeros(1, imgsrc.length());
Y = zeros(1, imgsrc.length());
R = zeros(1, imgsrc.length());
p = 256; % crop parameter
r1 = 20; % inner radius in pixels
r2 = 40; % outer radius in pixels
tic;
for t=1:imgsrc.length(); 
    im1 = double(imgsrc.imread(t,1));
    w = size(im1,1);
    cim1 = im1(w/2-p:w/2+p-1,w/2-p:w/2+p-1);    
    [x, y] = detectspots(cim1, 2, 1e-8);
    [v1,b1] = measure_intensity(cim1,x,y,r1,r2,0.1,3);
        
    im2 = double(imgsrc.imread(t,2));
    cim2 = im2(w/2-p:w/2+p-1,w/2-p:w/2+p-1); 
    [v2,b2] = measure_intensity(cim2,x,y,r1,r2,0.05,3);
        
    x = x + w/2 - p;
    y = y + w/2 - p;
    [X(t),Y(t)] = imgsrc.tostagecoords(x,y,t,1);      
    
    R(t) = (v1.MeanIntensity - b1.MeanIntensity) / (v2.MeanIntensity - b2.MeanIntensity);
    if mod(t,10) == 0
        showgraphs(imgsrc, im1, im2, t, X, Y,x,y,r1,r2,R);
    end    
    %mask = zeros(size(cim1));
    %mask(b1.PixelIdxList) = cim1(b1.PixelIdxList);    
    %imshow(mask,[]);
    %drawnow;
end
toc

end

function showgraphs(imgsrc, im1, im2, t, X, Y,x,y,r1,r2,R)
    subplot(121)
    ref1 = imgsrc.imref2d(t,1);
    ref2 = imgsrc.imref2d(t,2);
    imshowpair(im1,ref1,im2,ref2,'falsecolor','ColorChannels', [1 2 0]);
    hold on
    plot(imgsrc.stageposition(1:2:end,1), imgsrc.stageposition(1:2:end,2));
    plot(X(t),Y(t),'r+');
    plot(X(1:t),Y(1:t),'g-');
    ex1 = x+r1*cosd(0:360);
    ex2 = x+r2*cosd(0:360);
    ey1 = y+r1*sind(0:360);
    ey2 = y+r2*sind(0:360);
    [ex1, ey1] = imgsrc.tostagecoords(ex1,ey1,t,1);
    [ex2, ey2] = imgsrc.tostagecoords(ex2,ey2,t,1);
    plot(ex1,ey1,'m');
    plot(ex2,ey2,'m--');    
    hold off    
    axis([X(t)-500,X(t)+500,Y(t)-500,Y(t)+500]);
    title(sprintf('Frame:%d %.2f %.2f', t, x, y));  
    subplot(122)
    plot(R);    
    drawnow  
end

