function colorplot(x,y,v)
if 0
    surface([x;x],[y;y],[zeros(size(x));zeros(size(x))],[v;v],...
        'facecol','no',...
        'edgecol','interp',...
        'linew', 2);
else   
    cmap = colormap('jet');
    colormap('gray');    
    idx = ceil((v -min(v)) / (max(v)-min(v)) * 255)+1;
    min(idx)
    max(idx)
    for i=1:length(x)
        plot(x(i),y(i),'.','color',cmap(idx(i),:));
        hold on
    end
    colormap('gray');
    axis ij
end