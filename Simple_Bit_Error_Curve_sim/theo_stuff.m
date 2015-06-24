L = 2;                              % Distance between bits
M = 2^L;
k = log2(M); % bits per symbol
N = 1e2;                            % Number of symbols in simulation
Nsample = 1;                        % Oversampling rate
Es = 2;                             % Energy per symbol
Eb = Es / k;                        % Energy per bit
EbNo_dB(:,1) = linspace(0, 10, 11); % Energy per bit to noise power spectral density ratio (Eb/No)
SNR_dB(:,1) = EbNo_dB + 10*log10(k);
EbNo_lin(:,1) = 10.^(EbNo_dB / 10); % Eb/No values in linear scale

% theorySer_QPSK = ((M-1)/M)*erfc(sqrt((3/(M^2-1))*(k*EbNo_lin)));
% simSer_QPSK = erfc(sqrt(0.5*(10.^(Es_N0_dB/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(Es_N0_dB/10))))).^2;
theoryBer_QPSK = (1/k)*erfc(sqrt((k/2)*EbNo_lin));
theorySer_QPSK =1-(1- theoryBer_QPSK).^2;
close all
figure
semilogy(Es_N0_dB,theoryBer_QPSK,'b.-');
hold on
semilogy(Es_N0_dB,theorySer_QPSK,'mx-');

grid on
legend('BER-QPSK', 'SER-QPSK');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for QPSK(4-QAM)')