% For QPSK 2 system

clear all
L = 2;                              % Distance between bits
M = 2^L;
k = log2(M); % bits per symbol
N = 1e4;                            % Number of symbols in simulation
Nsample = 1;                        % Oversampling rate
Es = 2;                             % Energy per symbol
Eb = Es / k;                        % Energy per bit
EbNo_dB(:,1) = linspace(0, 10, 11); % Energy per bit to noise power spectral density ratio (Eb/No)
SNR_dB(:,1) = EbNo_dB + 10*log10(k);
EbNo_lin(:,1) = 10.^(EbNo_dB / 10); % Eb/No values in linear scale

% Define Arrays
% ------------
diffInxArray = zeros(length(EbNo_lin),N);
diffMinArray = zeros(length(EbNo_lin),N);
rxIdx = zeros(length(EbNo_lin),N);
rxSymbAssignment = zeros(length(EbNo_lin),N);
rxBits  = zeros(k,N);
bit_err = zeros(length(EbNo_lin),1);
diffCellArray = cell(length(EbNo_lin),1);
meanrxSymb = zeros(length(EbNo_lin),1);
meanAWGN = zeros(length(EbNo_lin),1);
meansymbArray = zeros(length(EbNo_lin),1);
meantxSymb = zeros(length(EbNo_lin),1);
meanNoise  = zeros(length(EbNo_lin),1);
diffAll = zeros(M,N);

% Symbol definition
% ------------
% symbArray =[-1-1i,1-1i,-1+1i,1+1i];
symbArray =[1,1i,-1,-1i];

% Define AWGN normalised
% ----------------------
AWGN = ((randn(1,N))+1i*(randn(1,N)))/sqrt(2);

% Random Data (binary to decimal)
% -------------------------------
txBits = reshape((rand(N*k,1) > 0.5),k,[]);                             % Random bit array for transmittion
txIdx = 2*txBits(1,:) + txBits(2,:);

% Modulation
% ------------
txSymb = symbArray(txIdx+1);     % From look up table determine symbol to be transmitted

for i = 1:length(EbNo_lin)
    % Adding AWGN to symbol
    rxSymb(i,:) = txSymb + sqrt(Eb/(EbNo_lin(i)))*AWGN;
    % Demodulation
    for ii = 1:N
        % Euclidean distance between recieved signal and symbol array
        diffAll(:,ii) = reshape(abs(symbArray-rxSymb(i,ii)).^2,length(symbArray),[]);
        [minDiff,diffIdx] = min(diffAll(:,ii));
        diffInxArray(i,ii) = diffIdx-1;
        rxSymbAssignment(i,ii) = symbArray(diffInxArray(i,ii)+1);     % From look up table determine symbol to be transmitted
        biStream = Bi2Dec( diffInxArray(i,ii), k );
        rxBits(:,ii) = reshape(biStream',k,[]);
    end
    bit_err(i) = sum(rxBits(k,:) ~= txBits(k,:));
end

% Bit Error
bit_err = bit_err / (k * N);
% bit_err_theo = (1/k)*erfc(sqrt(k*EbNo_lin));
bit_err_theo = (1/k)*erfc(sqrt((k/2)*EbNo_lin));
% Plot Results
hold on
figure(2)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theory', 'Simulation QPSK');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK Gray Code modulation');