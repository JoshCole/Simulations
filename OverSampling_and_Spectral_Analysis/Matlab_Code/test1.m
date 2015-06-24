clear all
dt = 0.001; n = 1000;
t = [0:n-1] * dt; f = 33;
s_t = cos( 2*pi*f*t );
S_w = abs(fftshift(fft(s_t)));
faxe = [-n/2:n/2-1]/n / dt;

figure
plot(t,s_t)
hold on
s_t = exp( 1i*2*pi*f*t );
plot(t,s_t,':r')
grid on
axis([0 0.2 -1 1])
hold off

figure
plot(faxe,S_w)
grid on
axis([-100 100 0 500])

% 
% S_w = abs(fftshift(fft(s_t)));
% 
% figure
% plot(t,s_t)
% grid on
% axis([0 0.2 -1 1])
% 
% figure
% plot(faxe,S_w)
% grid on
% axis([-100 100 0 500])