function colorplot(x,y,v,mode)
if mode == 1
    surface([x;x],[y;y],[zeros(size(x));zeros(size(x))],[v;v],...
        'facecol','no',...
        'edgecol','interp',...
        'linew', 2);
else   
    cmap = colormap('jet');
    colormap('gray');    
    idx = floor((v -min(v)) / (max(v)-min(v)) * (size(cmap,1)-1))+1;    
    for i=1:length(x)-1        
        line([x(i);x(i+1)],[y(i);y(i+1)],'color', cmap(idx(i),:));
        hold on
    end
    colormap('gray');
    axis ij
end