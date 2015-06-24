% Modulation Schemes
clear all
close all

modulationScheme = 'BPSK';
parity = 'gCRC8';
N = 1e6;                            % Number of symbols in simulation
Nsample = 1;                        % Oversampling rate
EbNo_dB = linspace(0, 10, 11); % Energy per bit to noise power spectral density ratio (Eb/No)
EbNo_lin = 10.^(EbNo_dB / 10); % Eb/No values in linear scale

switch modulationScheme
    case 'BPSK'
        L = 1;                              % Distance between bits
        M = 2^L;
        Es = 1;                             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[-1,1];
    case 'QPSK'
        L = 2;                              % Distance between bits
        M = 2^L;
        Es = 1;                             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[1,1i,-1,-1i];
    case 'QPSK-GrayCode'
        L = 2;                              % Distance between bits
        M = 2^L;
        Es = 2;                             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[-1-1i,1-1i,-1+1i,1+1i];
    case '16-QAM'
        A= 1;
        L = 4;                              % Distance between bits
        M = 2^L;
        Es = ((2*(M-1))/3)*A^2;             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[-3+3i,-3+1i,-3-3i,-3-1i,-1+3i,-1+1i,-1-3i,-1-1i,3+3i,3+1i,3-3i,3-1i,1+3i,1+1i,1-3i,1-1i];
end

k = log2(M);                        % bits per symbol
Eb = Es / k;                        % Energy per bit
SNR_dB = EbNo_dB + 10*log10(k);

% Transmit Data
% ------------
[ BER_10, BLER_10, CRCER_10 ] = transPacket( N, k, 10, Eb, EbNo_lin, symbArray, parity );
% [ BER_100, BLER_100, CRCER_100 ] = transPacket( N, k, 100, Eb, EbNo_lin,
% symbArray, parity );

% Bit Error Theory
% ---------------
switch modulationScheme
    case 'BPSK'
        bit_err_theo = (1/2)*erfc(sqrt(2*EbNo_lin)/sqrt(2));
    case 'QPSK'
        bit_err_theo = (1/k)*erfc(sqrt(k*EbNo_lin)/sqrt(2));
    case 'QPSK-GrayCode'
        bit_err_theo = (1/k)*erfc(sqrt(k*EbNo_lin)/sqrt(2));
    case '16-QAM'
        bit_err_theo = (1/k)*(3/2)*erfc(sqrt((k/10)*EbNo_lin));
end

figure(1)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,BER_10,'ko-');
grid on

legend('Theoretical bit error', 'Simulated bit error');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK modulation');