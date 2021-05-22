clc
clear all
close all

L =  512; %Length of filter
Fs = 6000; %Sampling frequency > 2*1633
fc = [697, 770, 852, 941, 1209, 1336, 1477, 1633]; 
len = 0:L-1;
t = 0:1/Fs:0.1-1/Fs;

%initialising zero matrices
h_n = zeros(8,L);
y_n = zeros(8,L+Fs/10-1);
H = zeros(8,L);
W = zeros(8,L);
pwr = zeros(1,8);

for ii = 1:8
    h_n(ii,:) = 0.0085*cos(len*(2*pi*fc(ii)/Fs)); %defining filter response in time domain
end

%Look-up table for symbols. Indices correspond to index of corresponding
%frequencies in fc array
%Arranging matrix like this according to MATLAB indexing convention
mat_symbol = ['1' '4' '7' '*'; '2' '5' '8' '0'; '3' '6' '9' '#'; 'A' 'B' 'C' 'D'];

symbol = input('Which symbol do you want to send?\n', 's');

%%% Encoder %%%
idx = find(mat_symbol == symbol);
%Access expressions to find corresponding frequencies
f1 = fc(floor((idx-1)/4)+1);
f2 = fc(4+mod(idx-1,4)+1);

% Generating the signal and adding noise
x = cos(2*pi*f1*t)+cos(2*pi*f2*t);
noise = randn(size(x)); %noise to simulate channel effects
noise = ((max(x)/max(noise))/10)*noise;
x = x + noise;
%%% Encoder End %%%

figure();
sgtitle('Input Signal in Time Domain')
subplot(211);
stem(t(1:50*floor(Fs/f2)),x(1:50*floor(Fs/f2)))
grid on
xlabel('Time in seconds');ylabel('x(t)');title('Input signal in time domain (discrete representation)');
subplot(212);
plot(t(1:50*floor(Fs/f2)),x(1:50*floor(Fs/f2)))
grid on
xlabel('Time in seconds');ylabel('x(t)');title('Input signal in time domain (continuous representation)');

figure();
freq = -Fs/2:Fs/L:Fs/2-Fs/L;
plot(freq, abs(fftshift(fft(x, L))))
grid on
xlabel('Frequency in Hz');ylabel('|X(\omega)|');title('Input Spectrum');

%%% Passing the input through filter bank %%% 
figure()
for ii = 1:8
    y_n(ii,:) = conv(x, h_n(ii,:));
    [H(ii,:), W(ii,:)] = freqz(y_n(ii,:), L);
    H(ii,:) = abs(H(ii,:));
    pwr(ii) = rms(y_n(ii,:))^2;
    subplot(4,2,ii)
    plot(freq, (abs(fftshift(fft(y_n(ii,:), L)))))
    grid on
    xlabel('Frequency in Hz');ylabel('|X(\omega)|');
    title("Output Spectrum for Bandpass filter centred at " + num2str(fc(ii)) + "Hz")
end
%%% End of filtering %%%

%%% Decoder %%%
[max, idx] = maxk(pwr, 2);
% Sorting to ensure that the higher and lower frequency components are in
% ascending order for access expressions to work correctly.
idx = sort(idx);
fprintf("The output power was highest for bandpass filters centered at %d Hz and %d Hz.\n", fc(idx(1)), fc(idx(2)));

% Access expressions based on chosen table
result = mat_symbol(idx(2)-4, idx(1));
fprintf("The chosen symbol, based on the output power of each filter, was %c.\n",result)

figure();
sgtitle('Output Signal in Time Domain')
subplot(211);
f1 = fc(idx(1)); f2 = fc(idx(2));
stem(t(1:100*floor(Fs/f2)),y_n(idx(1),1:100*floor(Fs/f2))+y_n(idx(2),1:100*floor(Fs/f2)))
grid on
xlabel('Time in seconds');ylabel('y(t)');title('Output signal in time domain (discrete representation)');
subplot(212);
plot(t(1:100*floor(Fs/f2)),y_n(idx(1),1:100*floor(Fs/f2))+y_n(idx(2),1:100*floor(Fs/f2)))
grid on
xlabel('Time in seconds');ylabel('x(t)');title('Output signal in time domain (continuous representation)');
%%% Decoder end %%%
