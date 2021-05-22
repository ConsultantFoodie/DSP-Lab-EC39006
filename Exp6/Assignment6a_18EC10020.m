clc
clear all
close all

order = 4; % order of all filters
[y,F] = audioread("./Samples/B1_A1.wav");
audiowrite("Original.wav",y, 10000);
[y,Fs] = audioread("Original.wav"); %reading the original speech signal
norm = Fs/2;

for N = [1,2,3,4] %Number of frequency bands
    for Fc = [16, 50, 160, 500] %Cutoff frequency (Hz) of envelope extraction filter
        [B_l, A_l] = butter(order, Fc/norm);
        noise = rand(size(y));
        output = zeros(size(y));
        envelope = zeros(size(y));

        %Band cut-off frequencies as decided by the paper
        bands = zeros(4,5);
        bands(1,:) = [1, norm-1, 0, 0, 0];
        bands(2,:) = [1, 1500, norm-1, 0, 0];
        bands(3,:) = [1, 800, 1500, norm-1, 0];
        bands(4,:) = [1, 800, 1500, 2500, norm-1];

        figure();
        sgtitle("Frequency Spectra of Noise Modulated Envelop (LPF Cut-off = "+Fc+" Hz)")
        for ii = 1:N
            [B, A] = butter(order, [bands(N,ii)/norm, bands(N,ii+1)/norm]);
            Y = filter(B,A,y); %Creating frequency bands
            Y_e = Y.*(Y>=0); %Half-wave rectification
            Y_el = filter(B_l, A_l, Y_e); %Envelope extraction by LPF
            noise_mod = noise.*Y_el; %Noise-modulation
            noise_mod = filter(B,A,noise_mod); %Spectral limitation

            subplot(N,1,ii);
            NUM = length(Y_el);
            f_range = -norm:2*norm/NUM:norm-1/NUM;
            plot(f_range, abs(fftshift(abs(fft(noise_mod))))); %Frequency spectrum of band
            xlabel("Frequency");ylabel("Magnitude");
            title("Frequency Spectrum of Band between "+bands(N,ii)+"Hz and "+bands(N,ii+1)+"Hz.");
            
            %Summing modulated noise and envelope for each band
            output = output + noise_mod; 
            envelope = envelope+Y_el;
        end
        %Passing output through LPF of cutoff 4kHz and amplifying
        [B_f, A_f] = butter(order, 4000/norm);
        output = filter(B_f, A_f, output);
        output = output.*(max(y)/max(envelope));
        %Passing envelope through LPF of cutoff 4kHz and amplifying
        envelope = filter(B_f, A_f, envelope);
        envelope = envelope.*(max(y)/max(envelope));
        t = 0:1/Fs:(length(y)-1)/Fs;
        
        %Plotting all signals in time domain
        figure();
        sgtitle("Signals and Envelopes for "+N+" bands and LPF cut-off = "+Fc+" Hz");
        subplot(311)
        plot(t, y);
        xlabel("Time");ylabel("Amplitude");
        title("Input Signal");

        subplot(312)
        plot(t, envelope);
        xlabel("Time");ylabel("Amplitude");
        title("Extracted Envelop");

        subplot(313);
        plot(t,output);
        xlabel("Time");ylabel("Amplitude");
        title('Envelop Modulated with Noise');

        %Writing output and envelope as audio files
        out_file = "./Audio_outputs/answer_"+N+"_freq_"+Fc+".wav";
        audiowrite(out_file,output,Fs);

        out_file = "./Audio_outputs/envelope_"+N+"_freq_"+Fc+".wav";
        audiowrite(out_file,envelope,Fs);
    end
end