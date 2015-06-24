close all
Fsub = 15e3;
Fs = 100*Fsub;
T = 1/Fsub;
t = linspace(0,T,300);
% t = t(1:end-1);
f1 = Fsub;
f2 = 2*Fsub;
phi = 2*pi*rand;
s1 = cos(2*pi*f1.*t+phi);
s2 = cos(2*pi*f2.*t+phi);
sumS = sum(s1.*s2);
figure
subplot(2,1,1)
plot(t,s1,'b')
hold on
plot(t,s2,'r')
xlabel('Time (sec)', 'FontSize',14);
ylabel('Amplitude', 'FontSize',14);
set(gca, 'FontSize',14); 
legend('15 kHz', '30 kHz')
axis([0 T -1 1])
subplot(2,1,2)
f = 0:Fsub:4.5e6-1;

stem(f,abs(fft(s1))/150)
hold on
stem(f,abs(fft(s2))/150,'r')
xlabel('Frequency (Hz)', 'FontSize',14);
ylabel('Amplitude', 'FontSize',14);
set(gca, 'FontSize',14); 
legend('15 kHz', '30 kHz')
axis([0 f2 0 1])