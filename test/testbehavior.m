%% test behavior identification
% load previously extracted tracks
clear all
addpath('../src')
tracks = trackset();


%% angle between two motion vector
smoothing = 2;
k = ceil(6 * smoothing + 1);
g = fspecial('gaussian',[1 k], smoothing);
shape = 'full';
[X,Y] = tracks.position();
d1p = conv(g, padarray([-1 1 0 0 0], [0 k-3]), shape);
v1 = [imfilter(X,d1p,'symmetric'); imfilter(Y,d1p,'symmetric')];
d1n = conv(g,  padarray([0 0 0 -1 1], [0 k-3]), shape);
v2 = [imfilter(X,d1n,'symmetric'); imfilter(Y,d1n,'symmetric')];
speed = sqrt(sum(v2.^2));
d2 = conv(g, padarray([1 -2 1], [0 k-3]), shape);
a = [imfilter(X,d2,'symmetric'); imfilter(Y,d2,'symmetric')];
% angle = acosd ( v1.v2 / |v1||v2| )
angle = acosd( sum(v1 .* v2) ./ (sqrt(sum(v1.^2)) .* sqrt(sum(v2.^2))));
% kappa =  |vx.ay-vy.ax|
kappa = abs(v1(1,:) .* a(2,:) - v1(2,:) .* a(1,:)) ./ max(sum(v1.^2),eps).^(3/2);
kappa(1) = kappa(2); kappa(end) = kappa(end-1);
figure(1)
clf, subplot(211), colorplot(X,Y,angle,1), axis equal
subplot(212), colorplot(X,Y,log(kappa),1), axis equal
%
%n=500
%plotyy(1:n,(log(kappa(1:n))),1:n,angle(1:n))

% relative velocity
v1 = v1 ./ repmat(sqrt(sum(v1.^2)),[2 1]);
v1n = [v1(2,:); -v1(1,:)];
vr = [dot(v2,v1); dot(v2,v1n)];
figure(2); clf; colorplot(vr(1,:),vr(2,:),R,1)
rho = sqrt(sum(vr.^2));
theta = atan2(vr(2,:),vr(1,:));
polar(theta,rho);
axis equal;box on; grid on
 %%
% get to direction from change of angle
angle_threshold = 90; % treshold
direction = mod(cumsum(double(angle>angle_threshold)),2);
% Exchange forward and backward if forward is less frequent.
if sum(direction==0) > sum(direction==1)
    direction = 1 - direction;
end
figure(3)
clf; colorplot(X,Y,direction,1);
box on; grid on;
xlabel('X [{\mu}m]');
ylabel('Y [{\mu}m]');

%%
subplot(211)
imagesc(direction);
subplot(212)
imagesc(speed)



