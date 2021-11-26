% define variables
fs = 1000
T = 1/fs
L = 1500
t = (0:L-1)*T

% perform fft
X = fft(x)
X = fftshift(X)

% double spectrum is desired
f = fs/L*(-L/2:L/2-1)

subplot(4, 1, 1)
plot(t, x, 'r')
xlabel('t [s]')
ylabel('x(t)')
title('Signal x(t)')

subplot(4, 1, 2)
plot(f, abs(X)/L)
title('Single-Sided Amplitude Spectrum of x(t)')
xlabel('f [Hz]')
ylabel('|X(f)|')

X(abs(X)/L < 0.13) = 0
subplot(4, 1, 3)
plot(f, abs(X)/L)
title('Filtered Single-Sided Amplitude Spectrum of x(t)')
xlabel('f [Hz]')
ylabel('|X(f)|')

x_filtered = ifft(ifftshift((X)), 'symmetric')
subplot(4, 1, 4)
plot(t, x_filtered)
xlabel('t [s]')
ylabel('x_{filtered}(t)')
title('Filtered Signal x(t)')




