% For QPSK system
clear all
close all

% System configuration
% ------------
b = 2;                              % Distance between bits
M = 2^b;
k = log2(M);                        % Bits per symbol
N = 1e4;                            % Number of symbols in simulatio
Es = 1;%2                             % Energy per symbol
Eb = Es / k;                        % Energy per bit
norValue = 2;
EbNo_dB(:,1) = linspace(0, 10, 11); % Energy per bit to noise power spectral density ratio (Eb/No)
EbNo_lin(:,1) = 10.^(EbNo_dB / 10); % Eb/No values in linear scale

% Define Arrays
% ------------
diffInxArray = zeros(length(EbNo_lin),N);
rxIdx = zeros(length(EbNo_lin),N);
bit_err = zeros(1,length(EbNo_lin));
diffCellArray = cell(length(EbNo_lin),1);
diffMinArray = zeros(length(EbNo_lin),N);
meanrxSymb = zeros(length(EbNo_lin),1);
meanAWGN = zeros(length(EbNo_lin),1);
meansymbArray = zeros(length(EbNo_lin),1);
meantxSymb = zeros(length(EbNo_lin),1);
meanNoise  = zeros(length(EbNo_lin),1);
rxBits  = zeros(k,N);

% Define AWGN
AWGN = ((randn(1,N))+1i*(randn(1,N)))/sqrt(2);

% Symbol definition
% ------------
% symbArray =[1,1i,-1,-1i];
symbArray =[-1-1i,1-1i,-1+1i,1+1i];
symbArrayNorm = mean(abs(symbArray).^2);
Esymb = symbArray./sqrt(symbArrayNorm);

%Modulation
% Random bit array for transmittion
txBits = reshape((rand(N*k,1) > 0.5),2,[]);

% Convert bits to Decimal
txIdx = 2*txBits(1,:) + txBits(2,:);

% Cordinates
% ----------
x = real(symbArray);
y = imag(symbArray);
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
        % Minimum Euclidean distance and associated index
        [diff,rxidx] = min(diffAll(:,ii));
        diffMin(i,ii) = diff;
        % Creates array of decimal index
        decIndex(ii)= rxidx-1;
        % Converts decimal into binary
        biStream = rem(floor(decIndex(ii)*pow2(1-k:0)),2);
        % Creates array of recieved bits
        rxBits(:,ii) = reshape(biStream',2,[]);
    end
    bit_err(i) = sum(rxBits(2,:) ~= txBits(2,:));
end

% Bit Error
% --------------------
bit_err = bit_err / (k * N);
bit_err_theo(:,1) = (1/k)*erfc(sqrt(EbNo_lin));

% Symbol Error
% --------------------
symb_err = sum(rxIdx(1,:) ~= txIdx(1,:));
symb_err = symb_err / (k * N);

% Check Results
% ------------
% for i = 1:length(EbNo_lin)
%     meantxSymb(i,1) = mean(abs(txSymb).^2);
%     meanrxSymb(i,1) = mean(abs(rxSymb(i,1:N)).^2);
%     meanAWGN(i,1) = mean(abs(AWGN(i,1:N)).^2);
%     meansymbArray(i,1) = mean(abs(symbArray).^2);
%     meanNoise(i,1) = mean(abs(Noise(i,1:N)).^2);
% end

% Plot Results
% ------------
figure(1)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theory', 'Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK modulation');