close all
clear all
Bandwith = 4.5e6;
deltaF = 15e3;
Nfft = 512;
N = 10;
T = 1/deltaF;
Ts = T/Nfft;
Nsubcarriers = 300;
TimeSample = (0:Nsubcarriers-1)*Ts;
Nwcdma = 256;
SubCarrierIndex_W_CDMA = [-Nwcdma/2:Nwcdma/2-1];

M = 4;
[ k, Es, Esnorm, Eb, Ebnorm, SymbArray ] = GetSymbolArrayData( M );
numBits = k*Nwcdma;

frame = [];

for numSymbols = 1:N
    Bits = rand(1,numBits) > 0.5;

    % Generate Symbol Index
    [ Idx, BlockSize ] = SymbolIndex( Bits, k );

    % Generate Symbols
    Symbols = SymbArray(Idx);
    
    % W-CDMA
    WDCMAymbol_QPSK = zeros(1,Nfft);


    % Map Symbols to Subcarrier
    WDCMAymbol_QPSK(SubCarrierIndex_W_CDMA+Nfft/2+1) = Symbols;

    frame = [frame WDCMAymbol_QPSK];
end

% Oversample
OversampleRate = 8;
[ InterperlationFilter ] = RaisedCosine( 0.5, 96 ,8 );
OVA = zeros(1,OversampleRate*length(frame));
OVA(1:OversampleRate:end) = frame;

x = conv(OVA,InterperlationFilter);

backoff = input('Enter backoff in dB > ');
ts = 1/(OversampleRate*7.68e6);
n = length(x);

GaindB = -15:1:5;
% GaindB = 0;
gn = length(GaindB);

reference_case_rms = 1.5237;   % Speech 12.2Kbps AMR cubic reference value
                                        % equivalent to 20*log10(v_rms) for
                                      % the Speech case
BW = [1.25, 2.5, 5.0, 10.0, 15.0, 20.0]*1e6;
perecentage = (60.*BW).*(0.015./BW);
BWocc = perecentage.*BW;
OBW = perecentage.*BW./1e6;
OBW_ref = 3.84;
OBPD_All = 10*log10(OBW/OBW_ref);
K = 1.56;
OBPD = OBPD_All(3);
for m = 1:gn
    xnew = zeros(1,n);
    y = zeros(1,n);
    for k=1:n
        G = 10^(GaindB(m)/10);
        xnew(k) = G.*x(k);
        [ AM_AM(k), AM_PM(k), y(k) ] = AmplifierModel( xnew(k),-1*backoff );
        xnew(k) = 10^(-1*backoff/20).*xnew(k);
    end
    
    v_abs=abs(y);
    v_scale=sqrt(mean(v_abs.*v_abs));
    v_normalised=v_abs/v_scale;

    v_cubed=v_normalised.*v_normalised.*v_normalised;
    v_rms=sqrt(mean(v_cubed.*v_cubed));
    raw_cubic_metric = 20*log10(v_rms);
    
    WPD = ((raw_cubic_metric-reference_case_rms)/K);
    cubic_OFDM = WPD+ OBPD;
    
    EVM = sqrt(mean((abs(y-xnew).^2))/mean(abs(xnew).^2));
    Peak_to_meanx = 10*log10(max(abs(x).^2)/mean(abs(x).^2));
    Peak_to_meany = 10*log10(max(abs(y).^2)/mean(abs(y).^2));
    disp(['Gain dB = ',num2str(GaindB(m))])
    disp(['EVM = ',num2str(EVM)])
    disp(['Peak to Mean dB x = ',num2str(Peak_to_meanx)])
    disp(['Peak to Mean dB y = ',num2str(Peak_to_meany)])
    disp(['Raw Cubic Metric dB  = ',num2str(raw_cubic_metric)])
    disp(['Cubic Metric dB = ',num2str(cubic_OFDM)])

    
    EVMArray(m) = EVM;
        
    N = length(xnew);
    Ex = sum(abs(xnew).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    xx(m) = XRMS;
    
    N = length(y);
    Ex = sum(abs(y).^2);
    Px = Ex/N;
    XRMS = sqrt(Px);
    yy(m) = XRMS;
    
    
    fftOFDM = 20*log10(abs(fft(y))*1/sqrt(length(y)));
    
    [ freq, psdy ] = WelchEstimate( y, Nfft, 0.5);
    welchy =10*log10(psdy);

    fs = 7.68*OversampleRate;
%     figure
%     subplot(2,1,1)
%     plot(0:fs/(length(x)-1):fs,fftOFDM); grid; title('Output of Amplifier');
%     ylabel('FFT in dB');
%     hold on
%     
% 
%     subplot(2,1,2)
%     
%     plot(0:fs/(length(welchy)-1):fs,welchy); grid; title('Output of Amplifier');
%     ylabel('Welch in dB'); xlabel('Frequency in Hz');
end

% [psdx,freq] = log_psd(x,n,ts);
% [psdy,freq] = log_psd(y,n,ts);
% [ freq, psdx ] = WelchEstimate( x, Nfft, ts );

% figure
pin = 10*log10(xx); % input power in dB
pout = 10*log10(yy); % output power in dB
% plot(pin,pin,'r'); 
% hold on
% plot(pin,pin-1,'g')
% plot(pin,pout); grid;
% xlabel('Input power - dB')
% ylabel('Output power - dB')
% legend('Input','1 dB Change', 'Output')
% 
% figure
% plot(pin,10*log10(abs(yy./xx)))
% grid on

figure
plot(pin,1-EVMArray)
grid on


