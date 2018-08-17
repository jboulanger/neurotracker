function y = tape(x, alpha)
% Tukey window

if nargin < 2
    alpha = .01;
end

if numel(size(x)) == 2
  wx = tape1D(size(x,1), alpha);
  wy = tape1D(size(x,2), alpha);
  w = repmat(wx', [1 size(x,2)]) .* repmat(wy, [size(x,1) 1]);
  y = w .* x;
else
    error('Not implemented')
end

function x = tape1D(N, alpha)
x = ones(1,N);
n = 0:N-1;
idx = find(n <= alpha * (N-1) / 2);
x(idx) = 0.5 * (1 + cos(pi *  (2 * n(idx) / (alpha * (N - 1) ) - 1)));
idx = find(n >= (N-1) * (1 - alpha/2));
x(idx) = 0.5 * (1 + cos(pi *  (2 * n(idx) / (alpha * (N - 1) ) - 2 / alpha + 1)));
