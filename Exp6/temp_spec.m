clc
clear all

%[y,F] = audioread("./Audio_Outputs/answer_4_freq_500.wav");
%[y1,F1] = audioread("./Audio_Outputs/answer_temp_4_freq_500.wav");
[y,F] = audioread("./Original.wav");
figure();
spectrogram(y);
%figure();
%spectrogram(y1);