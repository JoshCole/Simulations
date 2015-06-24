close all
clear all
% defining a signal in frequency domain
% subcarrier +1 alone
f = 7.68;
deltaf = 15e3;
T = 1/deltaf;
T_us = T*1e6;
Nfft = 512;
Nsubcarriers = 300;
t = 0:(T*1e6)/Nfft:T_us;
t = t(1:end-1);

OversampleRate = 1;
alpha = 0.22;                                                   % Roll Of Rate
truncation = 12;                                                % Truncation in respect to number of symbols
[ InterperlationFilter ] = RaisedCosine( alpha, truncation ,OversampleRate );      % Interperlation Filter


SubCarrierIndex = [-Nsubcarriers/2:Nsubcarriers/2-1];
% Real sine
figure
xtReal = sin(2.*pi.*1.*[0:1/Nfft:0.999]);
plot(t,xtReal)

meanSquareValue = xtReal*xtReal'/length(xtReal);
peakValue = max(xtReal.*conj(xtReal));
paprReal = peakValue./meanSquareValue;



% OVA_OFDMSymbol = zeros(OversampleRate*Nfft,1);
% OVA_OFDMSymbol(1:OversampleRate:end) = OFDMSymbol;
% 
% % Apply Interperation Filter
% Pulse_shapped = conv(OVA_OFDMSymbol,InterperlationFilter);
% 
% meanSquareValue = mean(abs(Pulse_shapped).^2);
% peakValue = max(abs(Pulse_shapped).^2);
% papr = peakValue./meanSquareValue;


% meanSquareValue = xt*xt'/length(xt);
% peakValue = max(xt.*conj(xt));
% papr = peakValue./meanSquareValue;
SubCarrier_1= find(SubCarrierIndex==1);
SubCarrier_2= find(SubCarrierIndex==2);
SubCarrier_3= find(SubCarrierIndex==3);
subcarriers = zeros(1,Nsubcarriers);


OFDMSymbol = zeros(1,Nfft);
% Map Symbols to Subcarrier

figure
subcarriers(SubCarrier_1)=1;
OFDMSymbol(SubCarrierIndex+Nfft/2+1) = subcarriers;
xt = ifft(fftshift(OFDMSymbol))*Nfft;
subplot(3,1,1)
plot(t,real(xt),'b')
hold on
plot(t,imag(xt),'g')
% plot(xt,'k','LineWidth',2)
axis([0 T_us min(imag(xt)) max(real(xt))])
set(gca, 'FontSize',14); 
xlabel('time {us}')
ylabel('amplitude')
title('complex sinusoidal')
legend('real', 'imag')
grid on
peakValue = max(xt.*conj(xt));
meanSquareValue = xt*xt'/length(xt);
papr = peakValue./meanSquareValue;

subplot(3,1,2)
subcarriers(SubCarrier_2)=1;
OFDMSymbol(SubCarrierIndex+Nfft/2+1) = subcarriers;
xt = ifft(fftshift(OFDMSymbol))*Nfft;
plot(t,real(xt),'b')
hold on
plot(t,imag(xt),'g')
axis([0 T_us min(imag(xt)) max(real(xt))])
set(gca, 'FontSize',14); 
xlabel('time {us}')
ylabel('amplitude')
title('complex sinusoidal')
legend('real', 'imag')
grid on
peakValue = max(xt.*conj(xt));
meanSquareValue = xt*xt'/length(xt);
papr = peakValue./meanSquareValue;

subplot(3,1,3)
subcarriers(SubCarrier_3)=1;
OFDMSymbol(SubCarrierIndex+Nfft/2+1) = subcarriers;
xt = ifft(fftshift(OFDMSymbol))*Nfft;
plot(t,real(xt),'b')
hold on
plot(t,imag(xt),'g')
axis([0 T_us min(imag(xt)) max(real(xt))])
peakValue = max(xt.*conj(xt));
meanSquareValue = xt*xt'/length(xt);
papr = peakValue./meanSquareValue;

set(gca, 'FontSize',14);  
xlabel('time {us}')
ylabel('amplitude')
title('complex sinusoidal')
legend('real', 'imag')
grid on
figure
freq = -f/2:f/(Nfft-1):f/2;
stem(freq,abs(fftshift(fft(xt,Nfft)/Nfft)))

% axis([0 512 0 1])
% 
%  fft_Symbols = fft(subcarriers,Nsubcarriers)/sqrt(Nsubcarriers);
%  SC_FDMASymbol = zeros(1,Nfft);
%  
% % Map Symbols to Subcarrier
% SC_FDMASymbol(SubCarrierIndex+Nfft/2+1) = fft_Symbols;
% % Apply IFFT
% SC_FDMA = ifft(SC_FDMASymbol,Nfft)*sqrt(Nfft);
% figure
% plot(t,real(SC_FDMA),'b')
% hold on
% plot(t,imag(SC_FDMA),'g')
% 
% meanSquareValue = SC_FDMA*SC_FDMA'/length(SC_FDMA);
% peakValue = max(SC_FDMA.*conj(SC_FDMA));
% papr = peakValue./meanSquareValue;