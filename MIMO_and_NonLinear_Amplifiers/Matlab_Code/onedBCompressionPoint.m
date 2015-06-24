close all
clear all
amplitudePointsdB = -15:1:3;
ampLin = 10.^(amplitudePointsdB / 10); 

output = zeros(1,length(ampLin));
for i = 1:length(ampLin)
    amp = ampLin(i);
    amp = repmat(amp,1,1000);
    % Amplifier Model
    [ AM_AM_OFDM_QPSK, AM_PM_OFDM_QPSK, xn ] = AmplifierModel( amp);
    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    output(i) = XRMS;
    
end
figure
% 
% hold on
plot(amplitudePointsdB,amplitudePointsdB,'g')
hold on
plot(amplitudePointsdB,amplitudePointsdB-1)
outdB =10*log10(abs(output));
plot(amplitudePointsdB,outdB,'r')
legend('Theoretical','1 dB Change', 'Amplifier Output')
figure
plot(amplitudePointsdB,10*log10(output./ampLin))
legend('Gain dB')
grid on

%%

amplitudePointsdB = -15:1:10;
ampLin = 10.^(amplitudePointsdB / 20); 

OnedBChange = -16:1:2;
output = zeros(1,length(ampLin));
input = zeros(1,length(ampLin));
Nfft = 128;
f1 = 8;

for i = 1:length(ampLin)
    amp = ampLin(i);
    MUx = ampLin(i).*exp(1i*2*pi*f1*[0:1/Nfft:0.999]);
    N = length(MUx);
    Ex = sum(abs(MUx).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    input(i) = XRMS;
    
    % Amplifier Model
    [ AM_AM_OFDM_QPSK, AM_PM_OFDM_QPSK, xn ] = AmplifierModel( MUx);
    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    output(i) = XRMS;
    
end
figure
% 
% hold on
plot(20*log10(abs(input)),20*log10(abs(input)),'g')
hold on
plot(20*log10(abs(input)),20*log10(abs(input))-1)
plot(20*log10(abs(input)),20*log10(abs(output)),'r')
legend('Theoretical','1 dB Change', 'Amplifier Output')

figure
plot(amplitudePointsdB,10*log10(output./ampLin))
legend('Gain dB')
%%
close all
clear all
amplitudePointsdB = -15:1:-3;


OnedBChange = -16:1:2;
output = zeros(1,length(OnedBChange));
input = zeros(1,length(OnedBChange));
Nfft = 128;
f1 = 8;
f2 = 13;
backoffdB = 0;

backoff = 10.^(backoffdB / 20); 
F = 0:Nfft-1;
input = zeros(1,length(amplitudePointsdB));
output = zeros(1,length(amplitudePointsdB));
for i = 1:length(amplitudePointsdB)
    
    ampLin = (10.^((amplitudePointsdB(i))/ 20)); 
    MUx = (ampLin.*exp(1i*2*pi*f1*[0:1/Nfft:0.999])+ampLin.*exp(1i*2*pi*f2*[0:1/Nfft:0.999]))/sqrt(2);
    MUx = backoff*MUx;
    
    N = length(MUx);
    Ex = sum(abs(MUx).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    input(i) = XRMS;
    
    % Amplifier Model
    [ AM_AM_OFDM_QPSK, AM_PM_OFDM_QPSK, xn ] = AmplifierModel( MUx);
    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    output(i) = XRMS;
    
    v = MUx;
    w = xn;
    e = abs(v-w).^2; 
    EVM(i) = sqrt(mean(e)/mean(abs(v).^2));   
%     figure
%     plot(real(MUx))
%     hold on
%     plot(real(xn),'r')
            
    MUxfft = 20*log10(abs((fft(MUx ,Nfft)*1/sqrt(Nfft))));
    CWfft = 20*log10(abs((fft(xn ,Nfft)*1/sqrt(Nfft))));
%     figure
%     plot(F, MUxfft)
%     hold on
%     plot(F,CWfft,'.-r');

    % Thrid Order 
    InputP0(i)=(MUxfft(f1+1)+MUxfft(f2+1))/2;
    f_IM3_1 = 2*f1-f2;
    f_IM3_2 = 2*f2-f1;
    P0(i) = (CWfft(f1+1)+CWfft(f2+1))/2;
    P03(i) = (CWfft(f_IM3_1+1)+CWfft(f_IM3_2+1))/2;
    
    IMD3 = P0(i)-P03(i);
    OIP3(i) = (IMD3/2)+P0(i);
    
    % 5th Order 
    f_IM5_2 = 3*f2-2*f1;
    P05(i) = CWfft(f_IM5_2+1);
    IMD5 = P0(i)-P05(i);
    OIP5(i) = (IMD5/4)+P0(i);
end
% Intercept Points
IIP3 = OIP3 - 10*log10(output./input);
IIP5 = OIP5 - 10*log10(output./input);

% Grad
x = InputP0;
y = IIP3;
m = zeros(1,length(OnedBChange)-1);
for k=0:length(y)-2
    m(k+1) = ((y(k+2) - y(k+1)))/(x(k+2) - x(k+1));
end

temp1 = InputP0(1);
temp2 = P03(1);
for i = 1:40
    P0new(i) = temp1;
    P03new(i) = temp2;
    temp1 = temp1+1;
    temp2 = temp2+3;
end

temp2 = P05(1);
for i = 1:40
    P05new(i) = temp2;
    temp2 = temp2+5;
end

figure
plot(20*log10(abs(input)),20*log10(abs(input)),'g')
hold on
plot(20*log10(abs(input)),20*log10(abs(input))-1)
plot(20*log10(abs(input)),20*log10(abs(output)),'r')
legend('Theoretical','1 dB Change', 'Amplifier Output')

% figure
% plot(P0new,P03new,'m')
% hold on
% plot(P0new,P0new,'g')

figure
semilogy(amplitudePointsdB,1-EVM)
xlabel('input dB')
ylabel('1-EVM')
grid on
input = (10.^((amplitudePointsdB)/ 20));
plot(amplitudePointsdB,10*log10(output./input))
legend('Gain dB')
grid on

%%
N = 10^6;

M = 64;
[ k, Es, Esnorm, Eb, Ebnorm, SymbArray ] = GetSymbolArrayData( M );
numBits = k*N;

Bits = rand(1,numBits) > 0.5;

% Generate Symbol Index
[ Idx, BlockSize ] = SymbolIndex( Bits, k );
        
% Generate Symbols
Symbols = SymbArray(Idx);

OversampleRate = 2;                                             % Oversample rate
alpha = 0.22;                                                   % Roll Of Rate
truncation = 32;                                                % Truncation in respect to number of symbols
[ InterperlationFilter ] = RaisedCosine( alpha, truncation ,OversampleRate );      % Interperlation Filter

% Oversample
OVA = zeros(1,OversampleRate*N);
OVA(1:OversampleRate:end) = Symbols;

% Apply Interperation Filter
Pulse_shapped = conv(OVA,InterperlationFilter);
        
    
backoffdB = -(1:10);
backoff = 10.^(backoffdB / 20); 
input = zeros(1,length(backoff));
output = zeros(1,length(backoff));
for i = 1:length(backoff)
    BF = backoff(i);
    xn = BF*Pulse_shapped;
    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    input(i) = XRMS;

    % Amplifier Model
    [ AM_AM_OFDM_QPSK, AM_PM_OFDM_QPSK, xn ] = AmplifierModel( xn);

    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    output(i) = XRMS;
end
inputdB = 20*log10(abs(input));
outputdB =20*log10(abs(output));

v_abs=abs(Pulse_shapped);
v_scale=sqrt(mean(v_abs.*v_abs));
v_normalised=v_abs/v_scale;

v_cubed=v_normalised.*v_normalised.*v_normalised;
v_rms=sqrt(mean(v_cubed.*v_cubed));
raw_cubic_metric = 20*log10(v_rms);

figure
plot(backoffdB,inputdB,'g')
hold on
plot(backoffdB,inputdB-1)
plot(backoffdB,outputdB,'r')
xlabel('Back Off dB')
ylabel('output dB')
legend('Theoretical','1 dB Change', 'Amplifier Output')
title(['Raw Cubic Metric', num2str(raw_cubic_metric)])

%%

N = 10^6;

M = 64;
[ k, Es, Esnorm, Eb, Ebnorm, SymbArray ] = GetSymbolArrayData( M );
numBits = k*N;

Bits = rand(1,numBits) > 0.5;

% Generate Symbol Index
[ Idx, BlockSize ] = SymbolIndex( Bits, k );
        
% Generate Symbols
Symbols = SymbArray(Idx);

OversampleRate = 2;                                             % Oversample rate
alpha = 0.22;                                                   % Roll Of Rate
truncation = 32;                                                % Truncation in respect to number of symbols
[ InterperlationFilter ] = RaisedCosine( alpha, truncation ,OversampleRate );      % Interperlation Filter

% Oversample
OVA = zeros(1,OversampleRate*N);
OVA(1:OversampleRate:end) = Symbols;

% Apply Interperation Filter
Pulse_shapped = conv(OVA,InterperlationFilter);
        
    
backoffdB = -(1:10);
backoff = 10.^(backoffdB / 20); 
input = zeros(1,length(backoff));
output = zeros(1,length(backoff));
for i = 1:length(backoff)
    BF = backoff(i);
    xn = BF*Pulse_shapped;
    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    input(i) = XRMS;

    % Amplifier Model
    [ AM_AM_OFDM_QPSK, AM_PM_OFDM_QPSK, xn ] = AmplifierModel( xn);

    N = length(xn);
    Ex = sum(abs(xn).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    output(i) = XRMS;
end
inputdB = 20*log10(abs(input));
outputdB =20*log10(abs(output));

v_abs=abs(Pulse_shapped);
v_scale=sqrt(mean(v_abs.*v_abs));
v_normalised=v_abs/v_scale;

v_cubed=v_normalised.*v_normalised.*v_normalised;
v_rms=sqrt(mean(v_cubed.*v_cubed));
raw_cubic_metric = 20*log10(v_rms);

figure
plot(backoffdB,inputdB,'g')
hold on
plot(backoffdB,inputdB-1)
plot(backoffdB,outputdB,'r')
xlabel('Back Off dB')
ylabel('output dB')
legend('Theoretical','1 dB Change', 'Amplifier Output')
title(['Raw Cubic Metric', num2str(raw_cubic_metric)])