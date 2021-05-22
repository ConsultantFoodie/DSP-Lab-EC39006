clc
clear all
close all

m = 8; %Filter order
mu = 1e-4; %Convergence rate (step size)
epsilon = 1e-3; % Error threshold

rng('default'); %For reproducible random numbers
cntr = 1; % For plotting all transfer functions in one figure

for k = [1,2,3,10] %%Choosing frequencies in kHz
f = k*1000;
fs = 4*f; %Sampling above Nyquist rate
t = 0:1/fs:0.1-1/fs;
N = length(t); % for N-point FFT
f_range = -fs/2:fs/N:fs/2-fs/N;

% Generating noisy sine wave 
x = 2*sin(2*pi*f*t)';
noise = randn(size(x));
x = x+noise;
x_n = buffer(x, m, m-1);
x_n = flip(x_n, 1);

% Initialising transfer function of ALE
h = zeros(m,1);
h_new = zeros(m,1);

% Passing the complete signal once
for ii = 1:size(x_n,2) %Passing vector of each sample
    h = h_new;
    y = x_n(:,ii)'*h; %Filtering (convolution)
    diff = x(ii)-y; %difference between desired and actual output
    update = x_n(:, ii)*diff;
    h_new = h + mu*update; %updating filter coefficients
    change = sum((h_new-h).^2)/sum(h.^2); %Finding the relative change
end
 
while change >= epsilon
    % Passing the complete signal till change >= error threshold
    for ii = 1:size(x_n,2) %Passing vector of each sample
        h = h_new;
        y = x_n(:,ii)'*h; %Filtering (convolution)
        diff = x(ii)-y; %difference between desired and actual output
        update = x_n(:, ii)*diff;
        h_new = h + mu*update; %updating filter coefficients
        change = sum((h_new-h).^2)/sum(h.^2); %Finding the relative change
    end
end

%Plotting transfer function of the filter
spectrum = fftshift(abs(fft(h, N)).^2);
figure(1);
sgtitle("Filter Order = "+m+", Step size (\mu) = "+mu+", Error (\epsilon) = "+epsilon);
subplot(2,2,cntr);
cntr = cntr+1;
plot(f_range, spectrum);
grid on;
xlabel('Frequency (Hz)'); title("Transfer function of adaptive filter for "+k+" KHz")
ylabel('$|H(e^{j2\pi f})|^2$', 'Interpreter', 'latex');

%Plotting input and output signals in time and frequency domain
figure();
sgtitle("Frequency of input sine wave = "+k+" kHz, Filter Order = "+m+", Step size (\mu) = "+mu+", Error (\epsilon) = "+epsilon);
subplot(2,2,1);
plot(t(1:20*fs/f),x(1:20*fs/f));
grid on;
xlabel("Time (seconds)"); ylabel("x(t)"); title("Input Signal in Time Domain");

subplot(2,2,2);
y = x_n'*h;
plot(t(1:20*fs/f),y(1:20*fs/f));
grid on;
xlabel("Time (seconds)"); ylabel("y(t)"); title("Output Signal in Time Domain");

subplot(2,2,3);
plot(f_range, fftshift(abs(fft(x))));
grid on;
xlabel("Frequency (Hz)"); ylabel("X(f)"); title("Input Signal in Frequency Domain");

subplot(2,2,4);
plot(f_range, fftshift(abs(fft(y))));
grid on;
xlabel("Frequency (Hz)"); ylabel("Y(f)"); title("Output Signal in Frequency Domain");
end