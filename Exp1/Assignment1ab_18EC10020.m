clc
clear all
close all

f = 1000; %frequency in Hz
fs = 12000; %sampling frequency in Hz
t = 0:1/fs:0.1; %time interval

x = 10*cos(2*pi*f*t)+6*cos(2*pi*2*f*t)+2*cos(2*pi*4*f*t); %Sampled Signal

subplot(3,2,1);
%Using plot to have an idea of what the continuous-time signal looks like
plot(t(1:6*fs/f), x(1:6*fs/f));
grid on;
title('Sampled x(t) as a continuous-time signal');
xlabel('Time(t) in seconds');
ylabel('|x(t)|');

subplot(3,2,2);
%Using stem rather than plot to show that this is a sampled signal
stem(t(1:6*fs/f), x(1:6*fs/f));
grid on;
title("x(t) sampled at " + fs/1000 + " KHz");
xlabel('Time(t) in seconds');
ylabel('|x(t)|');


plot_num = 3;
for N = [12, 64, 128, 256] %Iterating over desired values of N for N-point FFT
    dft = fftshift(fft(x, N)); % N-point DFT
    mag = abs(dft)/N; %Dividing by N to get similar energy distribution for all N
    f1 = -fs/2:fs/N:(fs/2 - fs/N); %frequency interval
    subplot(3,2,plot_num);
    plot_num = plot_num + 1;
    stem(f1, mag);
    grid on;
    title("DFT of sampled x(t) with N = " + N);
    xlabel('Frequency(f) in Hz');
    ylabel('|X(f)|/N');
end
