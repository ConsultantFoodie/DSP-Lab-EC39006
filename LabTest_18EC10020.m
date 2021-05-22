%%% Name: Hardik Suryaprakash Tibrewal
%%% Roll No.: 18EC10020
%%% DoB: 18/11/2000
clc
clear all
close all

%This is the recording of the audio part
%{ 
disp('Start speaking');
recObj = audiorecorder;
recordblocking(recObj, 7);
disp('End on recording.');

voice = getaudiodata(recObj);
audiowrite("LabAudio.wav", voice, 8000);
%}

%Reading the audio and plotting in time and frequency domains
[y, Fs] = audioread("LabAudio.wav");
figure();
t = 0:1/Fs:(numel(y)/Fs)-1/Fs;
plot(t,y);
grid on;
xlabel('Time (t) in seconds'); ylabel('Amplitude of y(t)');
title('Audio signal in time domain');

figure();
f_range = -Fs/2:Fs/numel(y):Fs/2-Fs/numel(y);
Y = fftshift(abs(fft(y)));
plot(f_range, Y);
grid on;
xlabel('Frequency (f) in Hz'); ylabel('Amplitude of Y(f)');
title('Audio signal in frequency domain');

%Making the required FIR filter using mentioned parameters
%Plotting frequency response of filter
M = mod((1+8+1+1+2+0+0+0),9);
f_cutoff = 150+10*M;
B = fir1(15, f_cutoff/(Fs/2), hamming(16));
figure();
sgtitle("Frequency and Phase response of the filter");
freqz(B,1,512);

%Filtering the signal and making plots
y_filt = filtfilt(B,1,y);
figure();
plot(t,y_filt);
grid on;
xlabel('Time (t) in seconds'); ylabel('Amplitude of Filtered y(t)');
title('Filtered Audio signal in time domain');

figure();
Y_filt = fftshift(abs(fft(y_filt)));
plot(f_range, Y_filt);
grid on;
xlabel('Frequency (f) in Hz'); ylabel('Amplitude of Filtered Y(f)');
title('Filtered Audio signal in frequency domain');

%Residual is the part that was removed by the filter. Finding and plotting
y_resid = y-y_filt;
figure();
plot(t,y_resid);
grid on;
xlabel('Time (t) in seconds'); ylabel('Amplitude of Filtered y(t)');
title('Residual Audio signal in time domain');

figure();
Y_resid = fftshift(abs(fft(y_resid)));
plot(f_range, Y_resid);
grid on;
xlabel('Frequency (f) in Hz'); ylabel('Amplitude of Residual Y(f)');
title('Residual Audio signal in frequency domain');

%Adding 25 dB average whie gaussian noise to the filter output. Plotting
%relevant plots
y_filt_noisy = awgn(y_filt, 25, 'measured', 'db');
figure();
plot(t,y_filt_noisy);
grid on;
xlabel('Time (t) in seconds'); ylabel('Amplitude of Noisy Filtered y(t)');
title('Noisy Filtered Audio signal in time domain');

figure();
Y_filt_noisy = fftshift(abs(fft(y_filt_noisy)));
plot(f_range, Y_filt_noisy);
grid on;
xlabel('Frequency (f) in Hz'); ylabel('Amplitude of Noisy Filtered Y(f)');
title('Noisy Filtered Audio signal in frequency domain');

%Ensuring the SNR is as desired
fprintf("SNR of the filtered signal is now = %.2f db\n", snr(y_filt, y_filt_noisy-y_filt));

%Writing audio file
audiowrite("18EC10020_OUTPUT.wav", y_filt_noisy, Fs);