close all
clear all
N1 = 512; Nsub1 = 300;
N2 = 2048; Nsub2 = 1200;
SubCarrierIndex1 = [-Nsub1/2:Nsub1/2-1];
SubCarrierIndex2 = [-Nsub2/2:Nsub2/2-1];

M = 4;
[ k_QPSK, Es, Esnorm, Eb, Ebnorm, SymbArray_QPSK ] = GetSymbolArrayData( M );
numBits_QPSK1 = k_QPSK*N1;
numBits_QPSK2 = k_QPSK*N2;

OversampleRate = 4;                                             % Oversample rate

alpha = 0.22;                                                   % Roll Of Rate
truncation = 12;                                                % Truncation in respect to number of symbols
[ InterperlationFilter ] = RaisedCosine( alpha, truncation ,OversampleRate );      % Interperlation Filter
NumberOfBlocks = 10^4;
Peak_to_mean_OFDM_QPSK1 = zeros(1,NumberOfBlocks);
Peak_to_mean_OFDM_QPSK2 = zeros(1,NumberOfBlocks);
raw_cubic_metric_QPSK1 = zeros(1,NumberOfBlocks);
raw_cubic_metric_QPSK2 = zeros(1,NumberOfBlocks);
for i = 1:NumberOfBlocks
        % Generate Bits
        Bits_QPSK1 = rand(1,numBits_QPSK1) > 0.5;
        % Generate Symbol Index
        [ Idx_QPSK1, BlockSize1 ] = SymbolIndex( Bits_QPSK1, k_QPSK );
        % Generate Symbols
        Symbols_QPSK1 = SymbArray_QPSK(Idx_QPSK1);
        
        % Generate Bits
        Bits_QPSK2 = rand(1,numBits_QPSK2) > 0.5;
        % Generate Symbol Index
        [ Idx_QPSK2, BlockSize2 ] = SymbolIndex( Bits_QPSK2, k_QPSK );
        % Generate Symbols
        Symbols_QPSK2 = SymbArray_QPSK(Idx_QPSK2);
        
        % OFDMA

        OFDMSymbol_QPSK1 = zeros(1,N1);
        OFDMSymbol_QPSK2 = zeros(1,N2);
        
        % Map Symbols to Subcarrier
        OFDMSymbol_QPSK1(SubCarrierIndex1+N1/2+1) = Symbols_QPSK1(1:Nsub1);
        OFDMSymbol_QPSK2(SubCarrierIndex2+N2/2+1) = Symbols_QPSK2(1:Nsub2);
        
        %Apply IFFT
        ifft_OFDMSymbol_QPSK1 = ifft(OFDMSymbol_QPSK1,N1)*sqrt(N1);
        ifft_OFDMSymbol_QPSK2 = ifft(OFDMSymbol_QPSK2,N2)*sqrt(N2);
        
        % Oversample
        ifft_OVA_QPSK1 = zeros(OversampleRate*N1,1);
        ifft_OVA_QPSK2 = zeros(OversampleRate*N2,1);
        
        ifft_OVA_QPSK1(1:OversampleRate:end) = ifft_OFDMSymbol_QPSK1;
        ifft_OVA_QPSK2(1:OversampleRate:end) = ifft_OFDMSymbol_QPSK2;
        
        % Apply Interperation Filter
        ifft_Pulse_shapped_QPSK1 = conv(ifft_OVA_QPSK1,InterperlationFilter);
        ifft_Pulse_shapped_QPSK2 = conv(ifft_OVA_QPSK2,InterperlationFilter);
        
        % Calculate Peak to Mean Ration
        Peak_to_mean_OFDM_QPSK1(1,i) = max(abs(ifft_Pulse_shapped_QPSK1).^2)/mean(abs(ifft_Pulse_shapped_QPSK1).^2);
        Peak_to_mean_OFDM_QPSK2(1,i) = max(abs(ifft_Pulse_shapped_QPSK2).^2)/mean(abs(ifft_Pulse_shapped_QPSK2).^2);
        
        filter_output = ifft_Pulse_shapped_QPSK1;
        v_abs=abs(filter_output);
        v_scale=sqrt(mean(v_abs.*v_abs));
        v_normalised=v_abs/v_scale;

        v_cubed=v_normalised.*v_normalised.*v_normalised;
        v_rms=sqrt(mean(v_cubed.*v_cubed));
        raw_cubic_metric_QPSK1(1,i) = 20*log10(v_rms);
        
        filter_output = ifft_Pulse_shapped_QPSK2;
        v_abs=abs(filter_output);
        v_scale=sqrt(mean(v_abs.*v_abs));
        v_normalised=v_abs/v_scale;

        v_cubed=v_normalised.*v_normalised.*v_normalised;
        v_rms=sqrt(mean(v_cubed.*v_cubed));
        raw_cubic_metric_QPSK2(1,i) = 20*log10(v_rms);
end
figure
plot(raw_cubic_metric_QPSK1)
hold on
plot(raw_cubic_metric_QPSK2,'r')
grid on

% amp = 1;
% F_s = 1000;
% F_w = 1;
% nsec = 2;
% dur= nsec*F_s;
% t = (0:1/F_s:nsec);
% n = 0:dur;
% ywave = amp*exp(1i*2*pi*n*F_w/F_s);
% 
% filter_output = ywave;
% v_abs=abs(filter_output);
% v_scale=sqrt(mean(v_abs.*v_abs));
% v_normalised=v_abs/v_scale;
% 
% v_cubed=v_normalised.*v_normalised.*v_normalised;
% v_rms=sqrt(mean(v_cubed.*v_cubed));
% raw_cubic_metric_QPSK1 = 20*log10(v_rms);
% 
% amp = 1;
% F_s = 1000;
% F_w = 10;
% nsec = 2;
% dur= nsec*F_s;
% t = (0:1/F_s:nsec);
% n = 0:dur;
% ywave = amp*exp(1i*2*pi*n*F_w/F_s);
% 
% filter_output = ywave;
% v_abs=abs(filter_output);
% v_scale=sqrt(mean(v_abs.*v_abs));
% v_normalised=v_abs/v_scale;
% 
% v_cubed=v_normalised.*v_normalised.*v_normalised;
% v_rms=sqrt(mean(v_cubed.*v_cubed));
% raw_cubic_metric_QPSK2 = 20*log10(v_rms);
% 
% figure
% plot(raw_cubic_metric_QPSK1)
% hold on
% plot(raw_cubic_metric_QPSK2,'r')
% grid on