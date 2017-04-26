function [region1,region2] = measure_intensity(img,x,y,r1,r2,p1,p2,smooth)
% [region1,region2] = measure_intensity(img,x,y,r1,r2,pvalue,smooth)
%
% Measure the intensities in two disks of radius r1 & r2 around position (x,y)
% Only a fraction of the pixel intensity is taken into account (pvalue)
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
filter = fspecial('gaussian', 3*smooth, smooth);
S = imfilter(double(img), filter, 'symmetric');

% distance function centered on (x,y)
[X,Y] = meshgrid(1:size(img,2),1:size(img,1));
d = sqrt((X-x).^2+(Y-y).^2);

% measure mean intensity for r < r1 & I > alpha
region1 = segment_bright_blob(S, d < r1, img, p1);
region2 = segment_background(S, d > r1 & d < r2, img, p2);

end

function region = segment_bright_blob(S, mask, img, pvalue)    
    mask = mask & (S > ecdfinv(S(mask), 1 - pvalue));        
    CC = bwconncomp(mask);    
    S = regionprops(CC,img, 'WeightedCentroid','PixelValues','PixelIdxList','MeanIntensity');
    I = arrayfun(@intensity, S);
    [~, idx] = max(I);
    region = S(idx);    
end

function region = segment_background(S,mask,img,pvalue)
    mask = mask & S > ecdfinv(S(mask), pvalue/2) & S < ecdfinv(S(mask), 1 - pvalue/2);    
    CC = bwconncomp(mask);
    S = regionprops(CC,img, 'WeightedCentroid','PixelValues','PixelIdxList','MeanIntensity');
    I = arrayfun(@intensity, S);
    [~, idx] = max(I);
    region = S(idx);    
end

function quantile = ecdfinv(X,pvalue)
    vals = sort(X);
    k = max(1,min(numel(vals), round(pvalue * numel(vals))));   
    quantile = vals(k);
end

function v = intensity(S)
    v = sum(S.PixelValues);
end

