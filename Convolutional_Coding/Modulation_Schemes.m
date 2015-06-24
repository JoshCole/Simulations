% Modulation Schemes
clear all
close all

modulationScheme = 'QPSK';
parity = 'gCRC24A';
% Generator Polynomials 3GPP TS25.212
% -----------------------------------
% genPoly1 = [1 1 1];
% genPoly2 = [1 0 1];
% 
% genPoly = [genPoly1;genPoly2];

genPoly1 = [1 0 1 1 0 1 1];
genPoly2 = [1 1 1 1 0 0 1];
genPoly3 = [1 1 1 0 1 0 1];

genPoly = [genPoly1;genPoly2;genPoly3];
memory = length(genPoly)-1;
% System configuration
% --------------------
N = 1e3;                            % Number of symbols in simulation
Nsample = 1;                        % Oversampling rate
EbNo_dB = 6;%linspace(0, 10, 11);      % Energy per bit to noise power spectral density ratio (Eb/No)
EbNo_lin = 10.^(EbNo_dB / 10);      % Eb/No values in linear scale

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

% Transmit and Recieve Data
% ------------
[ BER_10, BLER_10, CRCER_10 ] = tx_and_rx( N, k, 10, Eb, EbNo_lin, symbArray, parity, genPoly, memory );
% [ BER_100, BLER_100, CRCER_100 ] = tx_and_rx( N, k, 100, Eb, EbNo_lin,
% symbArray, parity, genPoly );

% Bit Error Theory
% ---------------
switch modulationScheme
    case 'BPSK'
        bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));
        Title = title('Block error probability curve for BPSK modulation');
    case 'QPSK'
        bit_err_theo = (1/k)*erfc(sqrt(k*EbNo_lin)/sqrt(2));
        Title = title('Block error probability curve for QPSK Gray Code modulation');
    case '16-QAM'
        bit_err_theo = (1/k)*(3/2)*erfc(sqrt((k/10)*EbNo_lin));
        Title = title('Block error probability curve for 16-QAM Gray Code modulation');
end

% Plot Results
% ---------------
figure(1)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,BER_10,'ko-');
% semilogy(EbNo_dB,BER_100,'co-');
semilogy(EbNo_dB,CRCER_10,'mx-');
% semilogy(EbNo_dB,CRCER_100,'rx-');
grid on

legend('Theory', 'Bit Error 10', 'Packet size 10');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');


figure(2)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,BLER_10,'ko-');
% semilogy(EbNo_dB,BLER_100,'co-');
semilogy(EbNo_dB,CRCER_10,'mx-');
% semilogy(EbNo_dB,CRCER_100,'rx-');
grid on

legend('Theory', 'Block Error 10', 'Packet size 10');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate')