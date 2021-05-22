clc
clear all
close all

N = 512; %declaring length of filter
k = floor((N-1)/2);
n = 0:1:(N-1);

wc = 0.8; %declaring normalised cut-off frequency in rads
w = -pi:1/2000:pi;

hd = zeros(1, N); %initialising ideal impulse response
for ii = 1:N
    if ii == k
        hd(ii) = wc/pi;
    else
        hd(ii) = sin(wc*(ii-k))/(pi*(ii-k));
    end
end

%declaring window functions
rectangular = ones(1, N);
triangular = 1 - 2*(n-(N-1)/2)/(N-1);
hanning = 0.5 - 0.5*cos((2*pi/(N-1))*n);
hamming = 0.54 - 0.46*cos((2*pi/(N-1))*n);
blackmann = 0.42 - 0.5*cos((2*pi/(N-1))*n) + 0.08*cos((4*pi/(N-1))*n);

%multiplication in time domain to get windowed filter
h_rect = hd.*rectangular;
h_trig = hd.*triangular;
h_hann = hd.*hanning;
h_hamm = hd.*hamming;
h_black = hd.*blackmann;

%%%% Generating an input signal and adding noise %%%%
f_pass = 500;
f_stop = 1500;
fs = 6000; %sampling frequency
t = 0:1/fs:(3*N-1)/fs;
rng('default') %to ensure same output from rand every time
noise = rand(1, 3*N);
x = 5*sin(2*pi*f_pass*t) + 5*sin(2*pi*f_stop*t);
add_noise = (max(x)/10)*noise/abs(max(noise)); %to control SNR
noisy_x = x + add_noise;
f_eq = -3000:2000/N:3000-2000/N;

h_matrix = [h_rect; h_trig; h_hann; h_hamm; h_black];

%%% Filtering the input signal through the various filters %%%
for ii=1:5
    if ii == 1
        name = "Rectangular";
    elseif ii == 2
        name = "Triangular";
    elseif ii == 3
        name = "Hanning";
    elseif ii == 4
        name = "Hamming";
    else
        name = "Blackmann";
    end
    y = filtfilt(h_matrix(ii,:), 1, x);
    y_n = filtfilt(h_matrix(ii,:), 1, noisy_x);
    
    %plotting time-domain response
    figure()
    sgtitle(name+" Window Time Domain");
    subplot(221);
    %plotting limited samples for clarity
    plot(t(1:floor(15*fs/f_stop)),x(1:floor(15*fs/f_stop)));
    grid on
    xlabel('Time in seconds'); ylabel('x(t)'); title('Input Signal');
    subplot(222);
    plot(t(1:floor(15*fs/f_stop)),noisy_x(1:floor(15*fs/f_stop)));
    grid on
    xlabel('Time in seconds'); ylabel('x(t)'); title('Input Signal with Noise');
    subplot(223);
    plot(t(1:floor(15*fs/f_stop)),y(1:floor(15*fs/f_stop)));
    grid on
    xlabel('Time in seconds'); ylabel('y(t)'); title('Output Signal');
    subplot(224);
    plot(t(1:floor(15*fs/f_stop)),y_n(1:floor(15*fs/f_stop)));
    grid on
    xlabel('Time in seconds'); ylabel('y(t)'); title('Output Signal with Noise');

    %plotting frequency domain response
    figure()
    sgtitle(name+" Window Frequency Domain");
    subplot(2,2,1)
    plot(f_eq, 20*log10(abs(fftshift(fft(x)))));
    xlabel('Frequency in Hz'); ylabel('Input Spectrum (in dB)'); title('Input Spectrum of Signal without Noise')
    grid on;
    subplot(2,2,2)
    plot(f_eq, 20*log10(abs(fftshift(fft(noisy_x)))));
    xlabel('Frequency in Hz'); ylabel('Input Spectrum (in dB)'); title('Input Spectrum of Signal with Noise')
    grid on;
    subplot(2,2,3)
    plot(f_eq, 20*log10(abs(fftshift(fft(y)))));
    xlabel('Frequency in Hz'); ylabel('Output Spectrum (in dB)'); title('Output Spectrum of Signal without Noise')
    grid on;
    subplot(2,2,4)
    plot(f_eq, 20*log10(abs(fftshift(fft(y_n)))));
    xlabel('Frequency in Hz'); ylabel('Output Spectrum (in dB)'); title('Output Spectrum of Signal with Noise')
    grid on;

    %printing amplitudes and SNR values
    max_out = max(y);
    max_out_noise = max(y_n-y);
    calc_snr = 20*log10(max_out/max_out_noise);
    fprintf("%0.2f %0.2f %0.2f \n",max_out, max_out_noise, calc_snr)
    disp("Input SNR(in dB) for "+name+" window = "+snr(x, noisy_x-x))
    disp("Output SNR(in dB) for "+name+" window = "+snr(y, y_n-y))
end
