function inspecttracks(obj)
%
% inspecttracks(obj)
%
% Inspect trackset object
%
% Jerome Boulanger 2017

f = figure('Visible','off',...
    'Name','Track image inspection tool','NumberTitle','off',...
    'Menubar', 'none', 'Toolbar', 'none');
N = size(obj.data,2);
sld = uicontrol('Style', 'slider',...
    'Min', 1, 'Max', N, 'Value', 1, 'SliderStep', [1/N, 100/N], ...
    'Position', [40 20 490 20], ...
    'Callback', @slider_callback, ...
    'Tag', 'slidertime', 'UserData', struct('tracks', obj));
txt = uicontrol('Style','text',...
        'Position',[200 45 120 20],...
        'String','Frame');
f.Visible = 'on';
end

function slider_callback(source, callbackdata)
h = findobj('Tag','slidertime');
data = h.UserData;
obj = data.tracks;
frame = int32(source.Value);
subplot(121)
obj.imshow(1, frame, 1);
hold on;
plot(obj.data(1,frame,1).x, obj.data(1,frame,1).y, 'w+','MarkerSize',10);
hold off
title(sprintf('ch1 frame:%d', frame))
subplot(122)
obj.imshow(1, frame, 2);
hold on;
plot(obj.data(1,frame,2).x, obj.data(1,frame,2).y, 'w+','MarkerSize',10);
hold off
title(sprintf('ch2 frame:%d', frame))
end