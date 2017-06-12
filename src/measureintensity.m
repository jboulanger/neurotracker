function [region1,region2] = measureintensity(img,x,y,radius,pvalue,smoothing)
% [region1,region2] = MEASUREINTENSITY(img,x,y,r1,r2,pvalue,smooth)
%
% Measure the intensities in two disks of radiuses [r1 r2] around position (x,y)
% Only a fraction of the pixel intensity is taken into account (pvalue = [p1,p2])
% 
% In the 1st ring r<r1, the average of 
%        intensity > cdfinv(p1) 
% is measured.
%
% In the 2nd ring r>r1 & r<r2, the average of
%    cdfinv(p2/2) < intensity < cdfinv(1-p2/2) 
% is measured.
% 
% Jerome Boulanger 2017

% pre-filter the image
if smoothing > 0
    filter = fspecial('gaussian', 3*smoothing, smoothing);
    S = imfilter(double(img), filter, 'symmetric');
else
    S = double(img);
end

% distance function centered on (x,y)
[X,Y] = meshgrid(1:size(img,2),1:size(img,1));
d = sqrt((X-x).^2+(Y-y).^2);

% measure mean intensity in inner disk
region1 = segment_bright_blob(S, d < radius(1), img, pvalue(1));
% measure mean intensity in outer ring
region2 = segment_background(S, d > radius(1) & d < radius(2), img, pvalue(2));
end

function region = segment_bright_blob(S, mask, img, pvalue)    
% Segment a bright blob in the mask based on the values of the smooth
% function. Return the intensity values in the original image img.        
    mask = mask & (S > prctile(S(mask), 100 * (1 - pvalue)));        
    CC = bwconncomp(mask);    
    allregion = regionprops(CC, img, 'WeightedCentroid', 'PixelValues','PixelIdxList','MeanIntensity','Area');
    I = [allregion.MeanIntensity];
    [~, idx] = max(I);
    region = allregion(idx);
    region.size = size(S);
    region.SMedianIntensity = median(S(region.PixelIdxList));
    region.MedianIntensity = median(img(region.PixelIdxList));
end

function region = segment_background(S,mask,img,pvalue)
% Segment the background in mask based on the values of the smooth
% function. Return the intensity values in the original image img.
    mask = mask & S > prctile(S(mask), 100*pvalue/2) & S < prctile(S(mask),100 * ( 1 - pvalue/2));    
    allregion = regionprops(mask, img, 'WeightedCentroid', 'PixelValues','PixelIdxList','MeanIntensity','Area');
    A = [allregion.Area];    
    [~, idx] = max(A);
    region = allregion(idx);
    region.size = size(S);
    region.SMedianIntensity = median(S(region.PixelIdxList));
end

function quantile = ecdfinv(X,pvalue)
% Empirical inverse cumultive distribution function
    vals = sort(X);
    k = max(1,min(numel(vals), round(pvalue * numel(vals))));   
    quantile = vals(k);
end


