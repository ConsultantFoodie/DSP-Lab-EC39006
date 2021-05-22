clc
clear all
close all

f = 1000; %frequency in Hz
fs = 20000; %sampling frequency in Hz
t = 0:1/(10*fs):0.1; %time interval to sample the signal at a much higher rate
ts = 0:1/fs:0.1; %time interval to sample to signal at desired rate of 20 kHz

x = square(2*pi*f*t); %A representation of the original signal
subplot(2,2,1);
plot(t(1:2000), x(1:2000));
grid on;
axis([0 0.01 -1.5 1.5]);
xlabel('Time in seconds'); ylabel('x(t)'); title('Original Signal x(t)');

xs = square(2*pi*f*ts); %Signal sampled at 20 kHz
subplot(2,2,2);
plot(ts(1:200), xs(1:200));
grid on;
axis([0 0.01 -1.5 1.5]);
xlabel('Time in seconds'); ylabel('x(t)'); title('x(t) sampled at 20 kHz');

plot_num = 3;
for N=[128, 256]
    Xs = fftshift(fft(xs, N)); %N-point DFT
    mag_s = abs(Xs)/N;
    f1 = -fs/2:fs/N:(fs/2 - fs/N); %frequency interval in Hz
    subplot(2,2,plot_num);
    plot_num = plot_num + 1;
    plot(f1, mag_s); %Using plot instead of stem for clarity
    grid on;
    xlabel('Frequency(f) in Hz'); ylabel('|X(f)|/N'); title("Fs = 20 kHz, N = " + N);
end
