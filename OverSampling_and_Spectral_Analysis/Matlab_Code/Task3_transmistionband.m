% Roll off effects
close all
clear all
Fs=8;
rollOfArray = [ 0 0.22  0.5 1];
n = length(rollOfArray);
H = zeros(193,n);
Nfft = 1024;
HF = zeros(Nfft,n);
for i = 1:n
    rollOff = rollOfArray(i);
    [ h ] = RootRaiseCosine( rollOff, 12, Fs );
    h = conv(h,h);
%     [ h ] = RaisedCosine( rollOff, 2^6, Fs );

    H(:,i) =h.';
    leg{i} = num2str(rollOff);
    
    Nh = length(h);
    zeroPad = zeros(Nfft -Nh,1);
    hf = [h.';zeroPad];
    HF(:,i)= fftshift(abs(fft(hf)));
end
n = length(h);
t = (-(n-1)/2:(n-1)/2)/Fs; 
figure
plot(t,H)
legend(leg)
xlabel('Time');
ylabel('Amplitude');
title('Impulse Response of Raised Cosine Filter');

figure
f = -Fs/2:Fs/(Nfft-1):Fs/2;
plot(f,20*log10(HF));
legend(leg)
xlabel('Frequency');
ylabel('Amplitude');
title('Frequency Response of Raised Cosine filter');