clc
clear all
close all

N = 8; %declaring window length
k = floor((N-1)/2);
n = 0:1:(N-1);

wc = 0.8; %normalized cut-off frequency in rads
w = -pi:1/2000:pi;

hd = zeros(1, N); %initialising desired impulse response
for ii = 1:N
    if ii == k
        hd(ii) = wc/pi;
    else
        hd(ii) = sin(wc*(ii-k))/(pi*(ii-k));
    end
end

%defining window functions
rectangular = ones(1, N);
triangular = 1 - 2*abs(n-(N-1)/2)/(N-1);
hanning = 0.5 - 0.5*cos((2*pi/(N-1))*n);
hamming = 0.54 - 0.46*cos((2*pi/(N-1))*n);
blackmann = 0.42 - 0.5*cos((2*pi/(N-1))*n) + 0.08*cos((4*pi/(N-1))*n);
impulse_mat = [hd; rectangular; triangular; hanning; hamming; blackmann];

%plotting ideal impulse and window functions
figure()
sgtitle("Impulse responses for N="+num2str(N));
for ii = 1:6
    if ii == 1
        name = "Ideal Impulse";
    elseif ii == 2
        name = "Rectangular Window";
    elseif ii == 3
        name = "Triangular Window";
    elseif ii == 4
        name = "Hanning Window";
    elseif ii == 5
        name = "Hamming Window";
    else
        name = "Blackmann Window";
    end
    subplot(3,2,ii);
    plot(n,impulse_mat(ii,:));
    xlim([0 N]);
    grid on
    xlabel('Time (or n)'); ylabel('x(n)'); title(name);
end

%multiplying windows with ideal impulse in time-domain
h_rect = hd.*rectangular;
h_trig = hd.*triangular;
h_hann = hd.*hanning;
h_hamm = hd.*hamming;
h_black = hd.*blackmann;

h_mat = [h_rect; h_trig; h_hann; h_hamm; h_black];

%plotting all frequency responses
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
    
    [H,W] = freqz(h_mat(ii,:), 1, w);
    figure()
    subplot(2,1,1)
    plot(W, abs(H));
    grid on
    xlabel('Normalized Frequency in rads/sample');
    ylabel('abs(H(w))');
    title("Frequency Response (Magnitude) for "+name+" Window when N=" + N);
    subplot(2,1,2)
    plot(W, angle(H));
    grid on
    xlabel('Normalized Frequency in rads');
    ylabel('angle(H(w)) in radians');
    title("Frequency Response (Phase) for "+name+" Window when N=" + N);    
end
