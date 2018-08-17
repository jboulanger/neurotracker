function y = ccwt(x,t,s)
% continuous wavelet transform with mexican hat
x = x(:);
t = t(:) - mean(t);
s = s(:);
y = zeros(numel(t),numel(s));
for j = 1:numel(s)
    g = 2 * (1 - (t/s(j)).^2) .* exp(-0.5 * (t /s(j)).^2) / (sqrt(3*s(j)) * pi^0.25);        
    y(:,j) = fftshift(ifft(fft(g).*fft(x)) / sqrt(s(j)));    
end