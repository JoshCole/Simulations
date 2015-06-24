close all
clear all
Fsub = 15e3;
T = 1/Fsub;
amp = sqrt(2/T);
F_s = 4.5e6;
F_w = 0:Fsub:4.5e6;
nsec = T;
dur= nsec*F_s;
t = (0:1/F_s:nsec);
n = 0:dur;
% ywave = amp*sin(2*pi*n*F_w/F_s);
F_w = F_w(1:5);
y = zeros(301,length(F_w));
ymax = 0;
for i = 1:length(F_w)
    phi = 0;%2*pi*rand;
    ywave = amp*exp(1i*2*pi*n*F_w(i)/F_s+phi);
    f = linspace(0,4.5e6,length(ywave));
    y(:,i) = ywave;
    if ymax> max(ywave)
        ymax = max(ywave);
    end
end
figure
subplot(2,1,1)
plot(t,y);
% title('Sine Wave 1 Hz, Two Period''s');
xlabel('Time (Sec)', 'FontSize',14);
ylabel('Amplitude', 'FontSize',14);
legend('0 kHz', '15 kHz', '30 kHz', '45 kHz','60 kHz')
set(gca, 'FontSize',14); 
subplot(2,1,2)
ffty = abs(fft(y))/301;
stem(f,ffty)
xlabel('Frequency (Hz)', 'FontSize',14);
ylabel('Amplitude', 'FontSize',14);
set(gca, 'FontSize',14); 
legend('0 kHz', '15 kHz', '30 kHz', '45 kHz','60 kHz')
axis([0 F_w(end) 0 max(ffty)])
%%
N= 300;
T= 1/15e3;
amp = sqrt(2/T);
Fsub = 4.5e6;
F = 4.5e6

Nsub = 300;
k = -Nsub/2:1:(Nsub/2);
k(k==0)=[];
% k = -Nsub/2:1:(Nsub/2)-1;
fs=2*Fsub;
nsec = 1/15e3;
dur = nsec*fs;
% t = 0:1/fs:nsec;
t = linspace(0,T,N);
G = zeros(N,Nsub);
for i = 1:Nsub
    phi = 0;%2*pi*rand;
    gkt = amp*exp((1i*2*pi*k(i).*t)/T+phi); 
    G(:,i)= gkt;
%     figure
%     plot(t,real(gkt))
end
figure
subplot(3,1,1)
fftG = abs(fft(G))/length(G);
% f = 0:Fsub/length(fftG)-1:4.5e6-1;
% f = 0:F/(length(fftG)-1):F;
%  f = linspace(-k(end)/T,k(end)/T,length(fftG));%-Fs/2:Fs/(n-1):Fs/2;
% F1 = -15e3*Nsub/2;
% F2 = 15e3*Nsub/2;
% F = 15e3*Nsub;
% f= F1:F/N:F2;
 f= -15e3*Nsub/2:15e3:15e3*Nsub/2;
 f(f==0) = [];

n = 6;
stem(f,fftG)
axis([-15e3*n 15e3*n 0 200])

subplot(3,1,2)
plot(t,real(G(:,1:1)))
subplot(3,1,3)
temp = unwrap(angle(G)*(180/pi));
plot(t,temp)

