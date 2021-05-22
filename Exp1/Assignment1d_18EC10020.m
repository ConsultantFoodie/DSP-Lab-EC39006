clc
clear all
close all

f = 1000; %frequency in Hz
fs1 = 12000; %initial sampling frequency
fs2 = 24000; %higher sampling frequency

t = 0:1/(10*fs1):0.1; %time interval to represent the original signal
ts1 = 0:1/fs1:0.1; %time interval to represent the signal sampled at 12 kHz
ts2 = 0:1/fs2:0.1; %time interval to represent the signal sampled at 24 kHz

x = 3*sin(2*pi*5*f*t)+5*sin(2*pi*3*f*t)+7*sin(2*pi*f*t); %Signal equation

figure();

len = (length(t)-1)/20;
plot(t(1:len), x(1:len)); %Only plotting first few samples for clarity
grid on;
axis([0, 0.005 -15 15]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('Original Signal (Sampled at 120 kHz)');

figure();

xs1 = 3*sin(2*pi*5*f*ts1)+5*sin(2*pi*3*f*ts1)+7*sin(2*pi*f*ts1); %12 KHz sampling

subplot(321);
len = (length(ts1)-1)/20;
stem(ts1(1:len), xs1(1:len)); %Only plotting first few samples for clarity
grid on;
axis([0, 0.005 -15 15]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('x(t) sampled at 12 kHz');

%Adding zeros between samples for interpolation
z = zeros(1,length(ts1));
xs2 = [xs1(:) z(:)]';
xs2 = xs2(:);

subplot(322);
len = (length(ts2)-1)/20;
stem(ts2(1:len), xs2(1:len)); %Only plotting first few samples for clarity
grid on;
axis([0, 0.005 -15 15]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('12 kHz samples with zero padding');

%Passing the samples through a 6th order butterworth filter of cut-off
%frequency of 6 kHz (normalized w.r.t half of sampling frequency
%as per function specification)
[b,a] = butter(6,(6*f)/(fs2/2));
filt_out = filter(b,a,xs2);

%plotting upsampled signal
up_xs1 = filt_out(1:length(filt_out)-1);
len = (length(ts2)-1)/20;
subplot(323);
stem(ts2(1:len), up_xs1(1:len)); %Only plotting first few samples for clarity
grid on;
axis([0, 0.005 -7.5 7.5]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('Upsampled x(t)');

len = (length(ts2)-1)/20;
subplot(325);
plot(ts2(1:len), up_xs1(1:len)); %Using plot rather than stem for clarity
grid on;
axis([0, 0.005 -7.5 7.5]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('Upsampled x(t)');

xs2 = 3*sin(2*pi*5*f*ts2)+5*sin(2*pi*3*f*ts2)+7*sin(2*pi*f*ts2); %24 kHz sampling

subplot(324);
len = (length(ts2)-1)/20;
stem(ts2(1:len), xs2(1:len)); %Only plotting first few samples for clarity
grid on;
axis([0, 0.005 -15 15]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('x(t) sampled at 24 kHz');

subplot(326);
len = (length(ts2)-1)/20;
plot(ts2(1:len), xs2(1:len)); %Using plot rather than stem for clarity
grid on;
axis([0, 0.005 -15 15]); %Adjusting axes
xlabel('Time(t) in seconds'); ylabel('x(t)'); title('x(t) sampled at 24 kHz');

%Comparing the upsampled signal with signal sampled at a higher rate in the
%frequency domain
figure();
N = 256;
up_ys1 = fftshift(fft(up_xs1, N)); %getting N-point DFT of upsampled signal
mag_ys1 = abs(up_ys1)/N;
f1 = -fs2/2:fs2/N:(fs2/2 - fs2/N); %frequency interval
subplot(1,2,1);
plot(f1, mag_ys1); %Using plot rather than stem for clarity
grid on;
xlabel('Frequency(f) in Hz'); ylabel('|X(f)|/N'); title('FFT of upsampled signal');

ys2 = fftshift(fft(xs2, N)); %getting N-point DFT of signal sampled at 24 kHz
mag_ys2 = abs(ys2)/N;
f1 = -fs2/2:fs2/N:(fs2/2 - fs2/N);
subplot(1,2,2);
plot(f1, mag_ys2); %Using plot rather than stem for clarity
grid on;
xlabel('Frequency(f) in Hz'); ylabel('|X(f)|/N'); title('FFT of signal sampled at 24 kHz');

