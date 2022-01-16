function y = ifft1d(x)
N = length(x);
z(2:N) = x(N:(-1):2);
z(1) = x(1);
y = (1/N)*fft1d(z);