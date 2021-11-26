% define vars
fs = 1000
T = 1/fs
L = 1500
t = (0:L-1)*T

% make figure
subplot(4, 1, 1)
plot(t, x, 'r')
xlabel('t [s]')
ylabel('x(t)')

% perform FFT
X = fft(x)
X = fftshift(X)

%X_shifted = fftshift(X)

% double spectrum is desired
f = fs/L*(-L/2:L/2-1)

% plot it
subplot(4, 1, 2)
plot(f, abs(X)/L)
title('Single-Sided Amplitude Spectrum of x(t)')
xlabel('f [Hz]')
ylabel('|X(f)|')

X(abs(X)/L < 0.13) = 0
subplot(4, 1, 3)
plot(f, abs(X)/L)
title('Single-Sided Amplitude Spectrum of x(t)')
xlabel('f [Hz]')
ylabel('|X(f)|')

x_filtered = ifft(ifftshift((X)), 'symmetric')
subplot(4, 1, 4)
plot(t, x_filtered)
xlabel('t [s]')
ylabel('x_{filtered}(t)')

% perform ifft
%x = ifft(P1)




