rect = ones(1,100);
No = 0.001;
awgn = sqrt(No)*randn(1,256)+1i*randn(1,256)/sqrt(2);
figure
plot(fftshift(abs(fft(rect,256))))

hold on
plot(abs(fft(awgn)),'r')