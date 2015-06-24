% For BPSK system
% -----------------
clear all
close all

% System configuration
% ------------
b = 1;                              % Distance between bits
M = 2^b;
k = log2(M);                        % Bits per symbol
N = 1e4;                            % Number of symbols in simulatio
A = 1;                              % Amplitude
Es = 1;                             % Energy per symbol
Eb = Es / k;                        % Energy per bit
EbNo_dB(:,1) = linspace(0, 10, 11); % Energy per bit to noise power spectral density ratio (Eb/No)
EbNo_lin(:,1) = 10.^(EbNo_dB / 10); % Eb/No values in linear scale

% Define Arrays
% ------------
diffInxArray = zeros(length(EbNo_lin),N);
rxIdx = zeros(length(EbNo_lin),N);
bit_err = zeros(1,length(EbNo_lin));
diffMinArray = zeros(length(EbNo_lin),N);
meanrxSymb = zeros(length(EbNo_lin),1);
meanAWGN = zeros(length(EbNo_lin),1);
meansymbArray = zeros(length(EbNo_lin),1);
meantxSymb = zeros(length(EbNo_lin),1);
meanNoise  = zeros(length(EbNo_lin),1);
rxBits  = zeros(k,N);
diffAll = zeros(M,N);

% Symbol definition
% ------------
symbArray =[-1,1];

% Define AWGN normalised
% ----------------------
AWGN = ((randn(length(EbNo_lin),N))+1i*(randn(length(EbNo_lin),N)))/sqrt(2);

% Random Data (binary to decimal)
% -------------------------------
txBits = reshape((rand(N*k,1) > 0.5),k,[]);     % Random bit array for transmittion
txIdx = txBits;                                 % Convert bits to Decimal

% Cordinates
% ----------
x = real(symbArray);
y = imag(symbArray);

% Modulation
% ------------
txSymb = symbArray(txIdx+1);     % From look up table determine symbol to be transmitted
Noise = repmat((sqrt(Eb./(EbNo_lin))),1,N) .* AWGN;

% Transmitt over Medium
% ------------
rxSymb = repmat(txSymb,length(EbNo_lin),1) + Noise;

for i = 1:length(EbNo_lin)
    for ii = 1:N
        % Euclidean distance between recieved signal and symbol array
        diffAll(:,ii) = reshape(abs(symbArray-rxSymb(i,ii)).^2,length(symbArray),[]);
        [minDiff,diffIdx] = min(diffAll(:,ii));
        diffInxArray(i,ii) = diffIdx-1;
        biStream = Bi2Dec( diffInxArray(i,ii), k );
        rxBits(:,ii) = reshape(biStream',k,[]);
    end
    bit_err(i) = sum(rxBits(:) ~= txBits(:));
end

% Bit Error
% --------------------
bit_err = bit_err / (k * N);
bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin));

% Check Results
% ------------
for i = 1:length(EbNo_lin)
    meantxSymb(i,1) = mean(abs(txSymb).^2);
    meanrxSymb(i,1) = mean(abs(rxSymb(i,1:N)).^2);
    meanAWGN(i,1) = mean(abs(AWGN(i,1:N)).^2);
    meansymbArray(i,1) = mean(abs(symbArray).^2);
    meanNoise(i,1) = mean(abs(Noise(i,1:N)).^2);
end

% Plot Results
% ------------
figure(1)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theoretical bit error', 'Simulated bit error');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK modulation');
