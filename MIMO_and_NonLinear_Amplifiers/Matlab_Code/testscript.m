clear
N = 10; % number of bits or symbols
Eb_N0_dB = [0:25]; % multiple Eb/N0 values
nTx = 2;
nRx = 2;
for ii = 1:length(Eb_N0_dB)

    % Transmitter
    ip = rand(1,N)>0.5; % generating 0,1 with equal probability
    s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0

    sMod = kron(s,ones(nRx,1)); % 
    sMod = reshape(sMod,[nRx,nTx,N/nTx]); % grouping in [nRx,nTx,N/NTx ] matrix

    h = 1/sqrt(2)*[randn(nRx,nTx,N/nTx) + j*randn(nRx,nTx,N/nTx)]; % Rayleigh channel
    n = 1/sqrt(2)*[randn(nRx,N/nTx) + j*randn(nRx,N/nTx)]; % white gaussian noise, 0dB variance

    % Channel and noise Noise addition
    y = squeeze(sum(h.*sMod,2)) + 10^(-Eb_N0_dB(ii)/20)*n;

    % Receiver

    % Forming the Zero Forcing equalization matrix W = inv(H^H*H)*H^H
    % H^H*H is of dimension [nTx x nTx]. In this case [2 x 2] 
    % Inverse of a [2x2] matrix [a b; c d] = 1/(ad-bc)[d -b;-c a]
    hCof = zeros(2,2,N/nTx)  ; 
    hCof(1,1,:) = sum(h(:,2,:).*conj(h(:,2,:)),1);  % d term
    hCof(2,2,:) = sum(h(:,1,:).*conj(h(:,1,:)),1);  % a term
    hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); % c term
    hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1); % b term
    hDen = ((hCof(1,1,:).*hCof(2,2,:)) - (hCof(1,2,:).*hCof(2,1,:))); % ad-bc term
    hDen = reshape(kron(reshape(hDen,1,N/nTx),ones(2,2)),2,2,N/nTx);  % formatting for division
    hInv = hCof./hDen; % inv(H^H*H)

    hMod =  reshape(conj(h),nRx,N); % H^H operation
    
    yMod = kron(y,ones(1,2)); % formatting the received symbol for equalization
    yMod = sum(hMod.*yMod,1); % H^H * y 
    yMod =  kron(reshape(yMod,2,N/nTx),ones(1,2)); % formatting
    yHat = sum(reshape(hInv,2,N).*yMod,1); % inv(H^H*H)*H^H*y
   
    % receiver - hard decision decoding
    ipHat = real(yHat)>0;

    % counting the errors
    nErr(ii) = size(find([ip- ipHat]),2);

end

simBer = nErr/N; % simulated ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 

close all
figure
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
semilogy(Eb_N0_dB,simBer,'mo-','LineWidth',2);
axis([0 25 10^-5 0.5])
grid on
legend('theory (nTx=1,nRx=1)', 'theory (nTx=1,nRx=2, MRC)', 'sim (nTx=2, nRx=2, ZF)');
xlabel('Average Eb/No,dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with 2x2 MIMO and ZF equalizer (Rayleigh channel)');