close all
clear all
% Create Reference Signal            
G = [1 0 0 0 0 1 1];
[ h ] = RootRaiseCosine( 0.22, 8 );
GF = G(2:end); 
N = length(G)-1; 
M = 2^N-1;
ProcessGainMsequence = 10*log10(M);
Reg = zeros(1,N); 
x = zeros(M,1);

% Initalise Reg 
Reg(end)=1;   
 
for n = 1:M
    buffer = mod(sum(GF.*Reg),2);      % Caluculate New Bit
    Reg = [buffer,Reg(1:N-1)];         % Shift Register   
    x(n) = Reg(N);
end
ReferenceSignal = (x.*2-1).';
MeanRef = mean(ReferenceSignal);
VarRef = var(ReferenceSignal);
CyclicPrefix = ReferenceSignal(M-14:M);
sig = [CyclicPrefix ReferenceSignal];
i = 1;
[ ReferenceSignalFreq, ReferenceSignalTime ] = CyclicCorr( sig(i:length(ReferenceSignal)),ReferenceSignal );
[Peak, Index] = max(ReferenceSignalTime);
figure
subplot(2,1,1)
stem(abs(ReferenceSignalTime))
title('Circular Correlation')
xlabel('Index');
ylabel('Magnitude');
i=Index;
test =sig(i:i+length(ReferenceSignal)-1);
[ ReferenceSignalFreq, ReferenceSignalTime ] = CyclicCorr(test,ReferenceSignal );
subplot(2,1,2)
stem(abs(ReferenceSignalTime))
title('Circular Correlation')
xlabel('Index');
ylabel('Magnitude');
%%
clear all
close all
% Create Reference Signal            
G = [1 0 0 0 0 1 1];
[ h ] = RootRaiseCosine( 0.22, 8 );
GF = G(2:end); 
N = length(G)-1; 
M = 2^N-1;
ProcessGainMsequence = 10*log10(M);
Reg = zeros(1,N); 
x = zeros(M,1);

% Initalise Reg 
Reg(end)=1;   
 
for n = 1:M
    buffer = mod(sum(GF.*Reg),2);      % Caluculate New Bit
    Reg = [buffer,Reg(1:N-1)];         % Shift Register   
    x(n) = Reg(N);
end
ReferenceSignal = (x.*2-1).';
MeanRef = mean(ReferenceSignal);
VarRef = var(ReferenceSignal);
CyclicPrefix = ReferenceSignal(M-14:M);
Nref = length(ReferenceSignal);
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( 4 );
TxPacket = double(rand(1,256) > 0.5);
Ntx = length(TxPacket);
[ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( TxPacket, k );

% Modulation, Transmit over Medium, Add Burst Error, Compute LLR
% --------------------------------------------------------------
[ TxSymbol ] = SymbolArray(TxSymbolIdx);    %  Gets Symbols to be Transmitted for Symbol Array
TxSymbolRe = reshape(TxSymbol,[],2);
sig_no_Delay = [TxSymbolRe(:,1).',CyclicPrefix, ReferenceSignal,TxSymbolRe(:,2).'];
% sig_no_Delay = [CyclicPrefix ReferenceSignal];
NyquistRate = 2;

% Up Sampled Signal 
% -----------------
[ TxUpSample, NTx] = OverSample( sig_no_Delay, NyquistRate );
[ TxRefOVS, NTxRef] = OverSample( ReferenceSignal, NyquistRate );
            
% Channel Delay
% -------------
SymbolyDelay = 15;
Delay = randi(SymbolyDelay,1)*NyquistRate; 
h = zeros(1,SymbolyDelay*NyquistRate);
h(Delay) = 1; % Delay Unit Power
SymbolDelay = (Delay/NyquistRate)-1; 
                     
% Convole transmitted signal will channel
sig = conv(h,TxUpSample);
 
Nrx = length(sig);
MidPoint = round(Nrx/2); 
iref = (80*NyquistRate);%MidPoint -round((Nref*NyquistRate)/2);
jref = iref+(Nref*NyquistRate)-1;
test = sig(iref:jref);
[ ReferenceSignalFreq, ReferenceSignalTime ] = CyclicCorr( test,TxRefOVS );
[Peak, Index] = max(ReferenceSignalTime);
figure
subplot(2,1,1)
stem(abs(ReferenceSignalTime))
 title(['Symbol Delay: ', num2str(SymbolDelay)])
xlabel('Index');
ylabel('Magnitude');
iref = iref + Index-1;
jref = iref+(Nref*NyquistRate)-1;
[ ReferenceSignalFreq, ReferenceSignalTime ] = CyclicCorr( sig(iref:jref),TxRefOVS );
subplot(2,1,2)
stem(abs(ReferenceSignalTime))
title('Circular Correlation')
xlabel('Index');
ylabel('Magnitude');