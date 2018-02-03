function x = psr_gauss_smoothing(x,dt,sigma)

% x: time-series
% dt: sampling period
% sigma: std of Gaussian

n = ceil(sigma / dt);
g = gausswin(n); % create Gaussian smoothing window
g = g / sum(g); % normalize
x = conv(x, g, 'same');

end