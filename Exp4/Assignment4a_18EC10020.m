clc
clear all
close all

mean = 0;
std_dev = 1;
N = 128;
rng('default');
noise = std_dev.*randn(N,1) + mean;
denom = 0.9;
A = [1, -denom, denom^2, -denom^3];

X = filter(1, A, noise); %Random Sequence
L = 8;

[H, W] = freqz(1,A,N/2);
H = (abs(H).^2).*std_dev^2;
H1 = [flip(H); H];
figure();
sgtitle("Transfer function: 1/(1 - "+abs(A(2))+"z^{-1} + "+abs(A(3))+"z^{-2} - "+abs(A(4))+"z^{-3})")
subplot(221);
f = -0.5:1/N:(0.5-(1/N));
plot(f, H1);
grid on;
xlabel('Frequency (f)');ylabel("Actual PSD");title("Known PSD using H(f)");

%%% Welch Method %%%
for ov_lap = [0, 0.25, 0.5] %fractional overlap
    M = floor(N/(L*(1-ov_lap)+ov_lap));
    D = floor(ov_lap*M);

    window = hamming(M);
    U = (1/M)*sum(window.^2);

    X_divs = zeros(M, L);
    P_xx = zeros(N, 1);
    l=1; r=M;

    for ii = 1:L
        X_divs(:,ii) = X(l:r);
        l = r-D+1;
        r = l+M-1;
    end

    for ii = 1:L
        X_divs(:,ii) = X_divs(:,ii).*window;
        P_xx(:) = P_xx(:) + (abs(fft(X_divs(:,ii),N)).^2)/(L*M*U);
    end

    subplot(2,2,ov_lap*4+2);
    plot(f, fftshift(P_xx));
    grid on;
    xlabel('Frequency (f)');ylabel("Estimated PSD");title("Estimated PSD using Welch Method "+num2str(100*ov_lap)+"% overlap");
end

figure();
sgtitle("Transfer function: 1/(1 - "+abs(A(2))+"z^{-1} + "+abs(A(3))+"z^{-2} - "+abs(A(4))+"z^{-3})")
subplot(221);
f = -0.5:1/N:(0.5-(1/N));
plot(f, H1);
grid on;
xlabel('Frequency (f)');ylabel("Actual PSD");title("Known PSD using H(f)");

%%% Parametric Method %%%
for p = [3, 6, 9]
    R = xcorr(X)/N; %Autocorrelation calculation
    r = R(N:N+p);
    
    %Matrix initialisation
    mat = zeros(p,p);
    mat2 = -1*(r(2:p+1));

    for ii = 1:p
        for jj = 1:p
            mat(ii,jj) = r(abs(ii-jj)+1);
        end
    end

    coeff_a = [1; mat\mat2]; %Coefficients. mat\mat2 gives inv(mat)*mat2
    new_var = sum(coeff_a.*r); %Variance

    [h_new, w_new] = freqz(1, coeff_a(:,1), N/2);
    h_new =  (abs(h_new).^2)*(new_var);
    l_new = [flip(h_new); h_new];
    subplot(2,2,p/3+1);
    plot(f,l_new);
    grid on;
    xlabel('Frequency (f)');ylabel("Estimated PSD");title("Estimated PSD using parametric method of order "+p);
end
