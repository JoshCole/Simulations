%%
% Simple Modulation - BPSK
% Clear all variables and Close all figures
clear all;                                            
close all;   

% For BPSK system
M = 2;
k = log2(M);

% Number of symbols in simulation
N = 1e6;

% energy per symbol
Es = 1;

% energy per bit (2 bits/symbol for QPSK)
Eb = Es / 2;

% Energy per bit to noise power spectral density ratio (Eb/No)
% Eb/No values to simulate at, in dB
% LINSPACE(X1, X2, N) generates N points between X1 and X2
EbNo_dB = linspace(0, 10, 11);

% Eb/No values in linear scale
EbNo_lin = 10.^(EbNo_dB / 10);

% Set array size for bit error
bit_err = zeros(size(EbNo_lin));

% Generating 0,1 with equal probability (randn(N,1) > 0.5)
Data = (randn(N,1) > 0.5);

% AWGN
temp = sqrt(Es/2)*(randn(1,N) + 1i*randn(1,N)); 
for ii = 1:N
    Noise(ii,1) = temp(1,ii);
end

for ii=1:length(EbNo_lin)
    % Generate source symbols, envelope as +/-1 
    Syms = (1 - 2 * Data);
    
    % Polar form
    r = abs(Syms);
    theta = angle(Syms)*180/pi;
    
    % Add Complex AWGN
    Syms_Noisy = Syms + sqrt(Eb/(2*EbNo_lin(ii))) * Noise;
    
    % recover symbols from each component (real and imaginary)
    syms_rec_r = sign(real(Syms_Noisy));
    
    % count bit errors
    bit_err(ii) = sum((syms_rec_r ~= real(Syms)));
end

% convert to bit error rate
bit_err = bit_err / (k*N);

% calculate theoretical bit error rate, functionally equivalent to:
% bit_err_theo = qfunc(sqrt(2*EbNo_lin));
bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));

% plot Results
close all
figure

semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theory', 'Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK modulation');
%%
% Simple Modulation - QPSK -1
% Clear all variables and Close all figures
clear all;                                            
close all;   

% For QPSK system
M = 4;
k = log2(M); 

% Number of symbols in simulation
N = 1e6;

% energy per symbol
Es = 1;

% energy per bit (2 bits/symbol for QPSK)
Eb = Es / 2;

% Energy per bit to noise power spectral density ratio (Eb/No)
% Eb/No values to simulate at, in dB
% LINSPACE(X1, X2, N) generates N points between X1 and X2
EbNo_dB = linspace(0, 10, 11);

% Eb/No values in linear scale
EbNo_lin = 10.^(EbNo_dB / 10);

% Set array size for bit error
bit_err = zeros(size(EbNo_lin));

for ii=1:length(EbNo_lin)
    
    % Generating 1,-1,1i,-1i with equal probability (randn(N,1) > 0.5)
    % Generate source symbols
    Syms =  1i * (1 - 2 * (randn(N, 1) > 0.5));
    
    % Polar form
    r = abs(Syms);
    theta = angle(Syms)*180/pi;
    
    % Add Complex AWGN
    Syms_Noisy = sqrt(Es/2) * Syms + sqrt(Eb/(2*EbNo_lin(ii))) * (randn(size(Syms)) + 1i * randn(size(Syms)));
    
    % recover symbols from each component (real and imaginary)
    syms_rec_r = sign(real(Syms_Noisy));
    syms_rec_i = sign(imag(Syms_Noisy));
    
    % count bit errors
    bit_err(ii) = sum((syms_rec_r ~= real(Syms)) + (syms_rec_i ~= imag(Syms)));
end
% convert to bit error rate
bit_err = bit_err / (k * N);

% calculate theoretical bit error rate, functionally equivalent to:
% bit_err_theo = qfunc(sqrt(2*EbNo_lin));
bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));

% plot Results
close all
figure

semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theory', 'Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK modulation');

%%
% Simple Modulation - QPSK -2
% Clear all variables and Close all figures
clear all;                                            
close all;   

% For QPSK system
M = 4;
k = log2(M); 

% Number of symbols in simulation
N = 1e6;

% energy per symbol
Es = 1;

% energy per bit (2 bits/symbol for QPSK)
Eb = Es / 2;

% Energy per bit to noise power spectral density ratio (Eb/No)
% Eb/No values to simulate at, in dB
% LINSPACE(X1, X2, N) generates N points between X1 and X2
EbNo_dB = linspace(0, 10, 11);

% Eb/No values in linear scale
EbNo_lin = 10.^(EbNo_dB / 10);

% Set array size for bit error
bit_err = zeros(size(EbNo_lin));

for ii=1:length(EbNo_lin)
    
    % Generating 0,1 with equal probability (randn(N,1) > 0.5)
    % Generate source symbols
    Syms = (1 - 2 * (randn(N,1) > 0.5)) + 1i * (1 - 2 * (randn(N, 1) > 0.5));
    
    % Add Complex AWGN
    Syms_Noisy = sqrt(Es/2) * Syms + sqrt(Eb/(2*EbNo_lin(ii))) * (randn(size(Syms)) + 1i * randn(size(Syms)));
    
    % recover symbols from each component (real and imaginary)
    syms_rec_r = sign(real(Syms_Noisy));
    syms_rec_i = sign(imag(Syms_Noisy));
    
    % count bit errors
    bit_err(ii) = sum((syms_rec_r ~= real(Syms)) + (syms_rec_i ~= imag(Syms)));
end
% convert to bit error rate
bit_err = bit_err / (k * N);

% calculate theoretical bit error rate, functionally equivalent to:
% bit_err_theo = qfunc(sqrt(2*EbNo_lin));
bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));

% plot Results
close all
figure

semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theory', 'Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK modulation');
%%
% Simple Modulation - QPSK -1
% Clear all variables and Close all figures
clear all;                                            
close all;   

% For QPSK system
M = 4;
k = log2(M); % bits per symbol

% Number of symbols in simulation
N = 1e6;

% energy per symbol
Es = 1;

% energy per bit (2 bits/symbol for QPSK)
Eb = Es / 2;

% Energy per bit to noise power spectral density ratio (Eb/No)
% Eb/No values to simulate at, in dB
% LINSPACE(X1, X2, N) generates N points between X1 and X2
EbNo_dB = linspace(0, 10, 11);

% Eb/No values in linear scale
EbNo_lin = 10.^(EbNo_dB / 10);

% Set array size for bit error
bit_err = zeros(size(EbNo_lin));

%symArray=[00,01,10,11]
symbArray =[1,1i,-1,-1i];

bitArray = reshape((randn(N*k,1) > 0.5),2,[]);

idx = 2*bitArray(1,:) + bitArray(2,:);
symb = symbArray(idx+1);
AWGN = sqrt(Eb/2)*((rand(1,N))+1i*(rand(1,N)));

for ii = 1:length(EbNo_lin)
    rxSymb = symb + EbNo_lin(ii)*AWGN;
end 

diffAll = (symbArray-rxSymb)^2;
[diff,rxidx] = min(diffAll);
