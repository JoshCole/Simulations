%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author        : Krishna Pillai
% Email         : krishna@dsplog.com
% Version       : 1.0
% Date          : 08 May 2009
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for computing the BER for BPSK modulation in Rayleigh
% channel with using rectangular pulse shaping at the transmitter
% and corresponding matched filter at the receiver

clear
N  = 10^6; % number of bits or symbols
T  = 1; % symbol duration of 1us 
os = 5; % oversampling factor
fs = 5/T; % sampling frequency in MHz

Eb_N0_dB = [0:10]; % multiple Eb/N0 values

for ii = 1:length(Eb_N0_dB)

   % Transmitter
   ip = rand(1,N)>0.5; % generating 0,1 with equal probability
   s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 

   % up sampling the signal for transmission
   sU = [s;zeros(os-1,length(s))];
   sU = sU(:).';
   sFilt = 1/sqrt(os)*conv(sU,ones(1,os));
   sFilt = sFilt(1:N*os);
   
   n = 1/sqrt(2)*[randn(1,N*os) + j*randn(1,N*os)]; % white gaussian noise, 0dB variance 
   
   % Noise addition
   y = sFilt + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

   % mathched filter
   yFilt = conv(y,ones(1,os)); % convolution
   ySamp = yFilt(os:os:N*os);  % sampling at time T
   
   % receiver - hard decision decoding
   ipHat = real(ySamp)>0;

   % counting the errors
   nErr(ii) = size(find([ip- ipHat]),2);

end

simBer = nErr/N; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
close all
figure
semilogy(Eb_N0_dB,theoryBer,'bs-','Linewidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-'),'Linewidth',2;
axis([0 10 10^-5 0.5])
grid on
legend('theory', 'sim-(matched filter)');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK modulation');

figure
[Pxx,W] = pwelch(sFilt,[],[],4096,fs,'twosided');
plot([-2048:2047]*fs/4096,10*log10(fftshift(Pxx1)),'b','LineWidth',2);
xlabel('frequency MHz');
ylabel('power spectral density');
title('Transmit spectrum - rectangular filtering, T = 1us')
grid on

