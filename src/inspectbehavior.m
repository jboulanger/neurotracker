function inspectbehavior(obj)
%
% inspecttracks(obj)
%
% Inspect trackset object
%
% Jerome Boulanger 2017

userdata = compute_userdata(obj);
scrsz = get(groot,'ScreenSize');
%f = figure('Visible','off',...
%    'Name','Behavior inspection tool','NumberTitle','off',...
%    'Menubar', 'none', 'Toolbar', 'none','Position',[1 1 scrsz(3) scrsz(4)]);

f = figure('Visible','off',...
    'Name','Behavior inspection tool','NumberTitle','off',...
    'Position',[1 1 scrsz(3) scrsz(4)]);

N = size(obj.data,2);
sld1 = uicontrol('Style', 'slider',...
    'Min', 1, 'Max', N, 'Value', 1, 'SliderStep', [1/N, 100/N], ...
    'Position', [100 30 scrsz(3)-200 20], ...
    'Callback', @slider1_callback, ...
    'Tag', 'slidertime', 'UserData', userdata);


txt = uicontrol('Style','text',...
        'Position',[100 10 scrsz(3)-200 20],...
        'String','Frame');
    
f.Visible = 'on';
refresh(userdata, 1);

end

function userdata = compute_userdata(obj)
[X,Y] = obj.position();
[R,T] = obj.ratio();
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
% normalize the speed by the time steps using the same derivative filter [um/s]
v2 = v2 ./ mean(diff(T)) * 1e-3; 
%W = log(abs(ccwt(diff(X,1)+1i*diff(Y,1),(0:numel(X)-2),logspace(0,3,64))))';
W = log(abs(ccwt(v1(1,:)+1i*v1(2,:),(0:numel(X)-1),logspace(0,3,64))))';
threshold = 50;
%score = (angle > threshold);
direction = mod(cumsum(double(angle>threshold)),2);
% Exchange forward and backward if forward is less frequent.
if sum(direction==0) > sum(direction==1)
    direction = 1 - direction;
end
% relative velocity (v2 proj on v1 and normal of v1)
v1n = [v1(2,:); -v1(1,:)];
userdata.relative_velocity = [dot(v2,v1); dot(v2,v1n)];
size(userdata.relative_velocity)
userdata.velocity = v2;
userdata.angle = angle;
userdata.curvature = curvature;
userdata.wavelet = W;
userdata.direction = direction;
userdata.R = R;
userdata.T = T;
userdata.X = X;
userdata.Y = Y;
end

function slider1_callback(source, callbackdata)
h = findobj('Tag','slidertime');
data = h.UserData;
%obj = data.tracks;
frame = int32(source.Value);
refresh(data, frame);
end

function refresh(obj, frame)

w = 1000;
t0 = max(1,frame-w);
t1 = min(numel(obj.angle),frame+w);

dirstr = {'backward','forward'};

subplot(241);
cla
colorplot(obj.X(t0:t1),obj.Y(t0:t1),obj.direction(t0:t1),1);
%plot(X(t0:t1),Y(t0:t1))
hold on;
plot(obj.X(frame), obj.Y(frame), 'ro', 'MarkerSize', 10);
quiver(obj.X(frame), obj.Y(frame), obj.velocity(1,frame),  obj.velocity(2,frame), 5);
hold off
grid on
box on
xlabel('X ({\mu}m)')
ylabel('Y ({\mu}m)')
title(sprintf('%s', dirstr{ obj.direction(frame)+1}));
axis square; axis ij

subplot(242)
plot(obj.angle(t0:t1)); hold on; plot([w+1 w+1],[0 max(obj.angle)],'r'); hold off
xlabel('time [s]');
title('angle [deg]')
axis square; axis tight; grid on;

subplot(243)
plot(obj.curvature(t0:t1));hold on; plot([w+1 w+1],[min(obj.curvature) max(obj.curvature)],'r'); hold off
title('curvature')
axis square; axis tight; grid on;

subplot(244)
speed = sqrt(sum(obj.velocity.^2));
plot(speed(t0:t1));hold on; plot([w+1 w+1],[min(speed) max(speed)],'r'); hold off
title('speed')
xlabel('time [s]');
ylabel('speed [um/s]');
axis square; axis tight; grid on;

subplot(245)
plot(obj.angle,speed,'-');
hold on;
plot(obj.angle(frame),speed(frame),'r.');
xlabel('angle [deg]');
ylabel('speed [um/s]');
hold off
axis square; axis tight; grid on;

subplot(246)
imagesc(obj.wavelet); hold on; plot([frame, frame],[0 64],'r');hold off
axis square; axis tight; grid on;
xlabel('time [s]');
ylabel('scale');
title('Wavelet transform');

subplot(247)
plot(obj.relative_velocity(1,:), obj.relative_velocity(2,:));
hold on;
plot(obj.relative_velocity(frame,1), obj.relative_velocity(frame,2),'ro');
hold off
axis square; axis tight; grid on;
xlabel('time [s]');
ylabel('scale');
title('Wavelet transform');

subplot(248)
plot(obj.R(t0:t1));hold on; plot([w+1 w+1],[min(obj.R),max(obj.R)],'r'); hold off
%plot(T,R,'r');%hold on; plot([w+1 w+1],[0,2],'r'); hold off
axis square; axis tight; grid on;
xlabel('time [s]');
ylabel('ratio');
title('Ratio');
end

