% For QPSK 2 system
clear all
close all

M = 16;
k = log2(M); % bits per symbol

% Number of symbols in simulation
N = 1e3;

% energy per symbol
Es = 1;

% Energy per bit
Eb = Es / k;

% Define AWGN
AWGN = ((randn(1,N))+1i*(randn(1,N)))/sqrt(2);

% Energy per bit to noise power spectral density ratio (Eb/No)
% Eb/No values to simulate at, in dB
% LINSPACE(X1, X2, N) generates N points between X1 and X2
EbNo_dB = linspace(0, 10, 11);

% Eb/No values in linear scale
EbNo_lin = 10.^(EbNo_dB / 10);

%Symbol definition
%symbArray =[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
%symbArray
%=[0000,0001,0011,0010,0110,0111,0101,0100,1100,1101,1111,1101,1010,1011,1001,1000];
symbArray =[1+1i,1+3i,3+1i,3+3i,1-3i,1-1i,3-3i,3-1i,-3+1i,-3+3i,-1+1i,-1+3i,-3-3i,-3-1i,-1-3i,-1-1i];

%Modulation
% Random bit array for transmittion
txBits = reshape((randn(N*k,1) > 0.5),k,[]);

% Convert bits to Decimal
idx = 8*txBits(1,:) + 4*txBits(2,:) + 2*txBits(3,:) + txBits(4,:);
% From look up table determine symbol to be transmitted
txSymb = symbArray(idx+1);

% Define Arrays
bit_err = zeros(size(EbNo_lin));
decIndex = zeros(1,N);
rxBits = zeros(2,N);
rxSymb = zeros(length(EbNo_lin),N);
diffAll = zeros(length(symbArray),N);
diffMin = zeros(length(EbNo_lin),N);


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
        biStream = rem(floor(decIndex(ii)*pow2(1-k:0)),k);
        % Creates array of recieved bits
        rxBits(:,ii) = reshape(biStream',k,[]);
    end
    bit_err(i) = sum(rxBits(k,:) ~= txBits(k,:));
end

% Bit Error
bit_err = bit_err / (k * N);
bit_err_theo = (1/k)*erfc(sqrt(k*EbNo_lin)/sqrt(2));

% Plot Results

figure(2)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,bit_err,'mx-');
grid on

legend('Theory', 'Simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for 16-QAM modulation');
