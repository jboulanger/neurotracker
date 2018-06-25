function inspectbehavior(obj)
%
% inspecttracks(obj)
%
% Inspect trackset object
%
% Jerome Boulanger 2017


f = figure('Visible','off',...
    'Name','Behavior inspection tool','NumberTitle','off',...
    'Menubar', 'none', 'Toolbar', 'none');
N = size(obj.data,2);
sld = uicontrol('Style', 'slider',...
    'Min', 1, 'Max', N, 'Value', 1, 'SliderStep', [1/N, 100/N], ...
    'Position', [5 5 480 20], ...
    'Callback', @slider_callback, ...
    'Tag', 'slidertime', 'UserData', struct('tracks', obj));
txt = uicontrol('Style','text',...
        'Position',[510 5 50 20],...
        'String','Frame');
f.Visible = 'on';
refresh(obj, 1);
end

function [v2, angle, curvature] = compute_descriptors(X,Y)
g = fspecial('gaussian',[1 11], 1);
shape = 'same';
d1p = conv(g, [0 0 0 0 -1 1 0 0 0 0 0], shape);
v1 = [imfilter(X,d1p,'symmetric'); imfilter(Y,d1p,'symmetric')];
d1n = conv(g, [0 0 0 0 0 -1 1 0 0 0 0], shape);
v2 = [imfilter(X,d1n,'symmetric'); imfilter(Y,d1n,'symmetric')];
d2 = conv(g, [0 0 0 0 1 -2 1 0 0 0 0], shape);
a = [imfilter(X,d2,'symmetric'); imfilter(Y,d2,'symmetric')];
% angle = acosd ( v1.v2 / |v1||v2| )
angle = acosd( sum(v1 .* v2) ./ (sqrt(sum(v1.^2)) .* sqrt(sum(v2.^2))));
% kappa =  |vx.ay-vy.ax|
curvature = abs(v2(1,:) .* a(2,:) - v2(2,:) .* a(1,:)) ./ max(sum(v2.^2),eps).^(3/2);
curvature(1) = curvature(2); curvature(end) = curvature(end-1);
curvature = log(curvature);
end

function slider_callback(source, callbackdata)
h = findobj('Tag','slidertime');
data = h.UserData;
obj = data.tracks;
frame = int32(source.Value);
refresh(obj, frame);

end

function refresh(obj, frame)
[X,Y] = obj.position();
w = 30;
[veloc, angle, curvature] = compute_descriptors(X,Y);
t0=max(1,frame-w);
t1=min(numel(angle),frame+w);
threshold = 50;
%score = (angle > threshold);
direction = mod(cumsum(double(angle>threshold)),2);
dirstr = {'backward','forward'};

subplot(241);
cla
colorplot(X(t0:t1),Y(t0:t1),direction(t0:t1),1);
%plot(X(t0:t1),Y(t0:t1))
hold on;
plot(X(frame), Y(frame), 'ro', 'MarkerSize', 10);
quiver(X(frame), Y(frame), veloc(1,frame), veloc(2,frame), 5);
hold off
grid on
box on
xlabel('X ({\mu}m)')
ylabel('Y ({\mu}m)')
title(sprintf('%s', dirstr{direction(frame)+1}));
axis square; axis ij

subplot(242)
plot(angle(t0:t1)); hold on; plot([w+1 w+1],[0 max(angle)],'r'); hold off
title('angle')
axis square; axis tight; grid on;

subplot(243)
plot(curvature(t0:t1));hold on; plot([w+1 w+1],[min(curvature) max(curvature)],'r'); hold off
title('curvature')
axis square; axis tight; grid on;

subplot(244)
speed = sqrt(sum(veloc.^2));
plot(speed(t0:t1));hold on; plot([w+1 w+1],[min(speed) max(speed)],'r'); hold off
title('speed')
axis square; axis tight; grid on;

subplot(245)
plot(angle,speed,'-');
hold on;
plot(angle(frame),speed(frame),'r.');
xlabel('angle [deg]');
ylabel('speed [um/s]');
hold off
axis square; axis tight; grid on;

subplot(247)
ratio = obj.ratio();
plot(ratio(t0:t1),'r');hold on; plot([w+1 w+1],[min(ratio) max(ratio)],'r'); hold off

axis square; axis tight; grid on;
title('Ratio');
end

