%HELP one_shot_ML_equalizer
%
%Illustrate how to use the ML space time code equalizer. This function simply
%simulate a MIMO communication with STBC encoding. Symbols are generated randomly from a specified constellation.
%Then, the receiver performs a Maximum Likelihood equalization to estimate the transmitted symbols
%
%Programmed by V. Choqueuse (vincent.choqueuse@gmail.com)

%Simulation parameters

N=512;                  %Number of symbols to be transmitted 
code_name='Alamouti';   %Space time code (see file space_time_coding to obtain the list of supported STBC)
rate='1';               %Space time code (see file space_time_coding to obtain the list of supported STBC)
num_code=1;             %Space time code (see file space_time_coding to obtain the list of supported STBC)
modulation='PSK';       %supported modulation PSK, QAM
state_nb=4;             %modulation with 4 states (4-PSK -> QPSK)
nb_receivers=4;         %Number of 4 receivers
snr=10;                  %Signal to noise ratio (dB)

fprintf('\n+-----------------------------------+\n');
fprintf('|   Simulate a MIMO communication   |\n');
fprintf('+-----------------------------------+\n\n');

%% extract space time block coding information
code_rate=str2num(rate);
[nb_emitters,code_length]=size(space_time_coding(0,code_name,rate,num_code,1));
Nb_symbole_code=code_length*str2num(rate);

%% Generate a symbol sequence randomly. The symbols belong to the set of
%%integer: [0 state_nb-1]
fprintf('- Generate %d random symbols: ',N);
M=4;
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
txBits = rand(1,code_rate*N*k) > 0.5;
[ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( txBits, k );
fprintf('\t\t\tOK\n');
  
%% Modulate the symbols
fprintf('- Apply %d-%s constellation: ',state_nb,M);    
[ modulated_symbols ] = SymbolArray(TxSymbolIdx);
fprintf('\t\t\tOK\n');

%% perform space time encoding
fprintf('- Perform %s-%s STBC encoding:',rate,code_name);
[STBC_blocs]=space_time_coding(modulated_symbols,code_name,rate,num_code);
fprintf('\t\tOK\n');

%% Create a random channel matrix
fprintf('- Generate a %d * %d Random Channel: ',nb_receivers,nb_emitters);
channel_matrix=sqrt(0.5)*(randn(nb_receivers,nb_emitters)+i*randn(nb_receivers,nb_emitters));
received_signal=channel_matrix*STBC_blocs;
fprintf('\t\tOK\n');

%% Apply AWGN noise
fprintf('- Apply %d dB additive noise: ',snr);
noise_variance=1/(10^(snr/10));
bruit=(sqrt(noise_variance/2))*(randn(nb_receivers,size(STBC_blocs,2))+...
                                i*randn(nb_receivers,size(STBC_blocs,2)));                          
received_signal=received_signal+bruit;
fprintf('\t\t\tOK (noise variance=%f)\n',noise_variance);

%% Perform Space time block equalization
fprintf('- Perform ML equalization: ');
modulator = SymbolArray;
equalized_symbols=coherent_ML_receiver(received_signal,channel_matrix,code_name,rate,num_code,modulator);
equalized_symbols=equalized_symbols(:).'; %convert matrix -> row vector
fprintf('\t\t\t\tOK\n');

%% Compare the equalized symbols with the real transmitted ones.
fprintf('- Compute Error rate: ');
LLRData  = LLR( LLRBitArray, equalized_symbols, k, noise_variance );
[RxSoftDecision, RxHardDecision] = DecisionType( LLRData, k, N );
estimated_symbols=demodulate(demodulator,equalized_symbols);    %demodulation
[snum_CSI,srate_CSI] = symerr(estimated_symbols,symbols);       %symbol error rate
[bnum_CSI,brate_CSI] = biterr(estimated_symbols,symbols);       %bit error rate
fprintf('\t\t\t\t\tSER= %f (%i error)\n\t\t\t\t\t\t\t\t\t\tBER= %f (%i error)\n',srate_CSI,snum_CSI,brate_CSI,bnum_CSI);