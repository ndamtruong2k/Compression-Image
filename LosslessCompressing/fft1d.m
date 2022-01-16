function y = fft1d(x)
N = length(x);
if N == 1
  y = x(1);
else
even = x(1:2:(N-1));
odd = x(2:2:N);
ye = fft1d(even);
yo = fft1d(odd);
D=exp(-2*pi*1j*(0:(N/2-1))'/N);
y = [ ye + yo.*D; ye - yo.*D];
end
