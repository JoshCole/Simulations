Fs = 8;
numBits = 2^16;
t  = (0:1/Fs:numBits-1/Fs).'; 
ovTxbit = zeros(Fs*numBits,1);
txbit = (rand(numBits,1)<0.5)*2-1;
ovTxbit(1:Fs:end) = txbit;
x= ovTxbit;
Nx = length(x);
Nfft = 512;
[ w ] = hanning( Nfft );
[Pxx,F] = pwelch(x,w,[],Nfft,Fs);
plot(F,Pxx)
%%
Fs = 1000;
t = 0:1/Fs:1-(1/Fs);
% 200Hz cosine + noise
x = cos(2*pi*t*200) + randn(size(t));
[Pxx,F] = pwelch(x,128,120,length(x),Fs,'onesided');
plot(F,10*log10(Pxx))