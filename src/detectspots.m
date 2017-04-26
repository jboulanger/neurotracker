function [x,y,v,pixels] = detectspots(img, scale, pvalue)
% return the coordinates and intensity of the brightest spots in an image
    filter = fspecial('gaussian', 3*scale, scale) - fspecial('gaussian', 3*scale, 3*scale);
    S = imfilter(double(img), filter, 'symmetric');
    m = median(S(:));
    s = median(abs(S(:)-median(S(:)))) * 1.4816;
    CC = bwconncomp(S > norminv(1 - pvalue, m, s));
    S = regionprops(CC,img, 'WeightedCentroid','PixelValues','PixelIdxList');
    I = arrayfun(@intensity, S);
    [v, idx] = max(I);   
    x = S(idx).WeightedCentroid(1);
    y = S(idx).WeightedCentroid(2);
    pixels = S(idx).PixelIdxList;
end

function v = intensity(S)
    v = mean(S.PixelValues);
end