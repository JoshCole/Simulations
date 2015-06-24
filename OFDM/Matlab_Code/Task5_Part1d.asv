% Cyclic Prefix
close all
amp = 1;
F_s = 4.5e6;
F_w = 15e3;
nsec = (1/F_w)*1;
dur= nsec*F_s;
t0 = 0:1/F_s:(35/F_s);
t = (t0(end):1/F_s:nsec+t0(end));
n = 0:dur;
y = amp*sin(2*pi*n*F_w/F_s);
CP =y(end-35:end); 
% y = amp*exp(1i*2*pi*n*F_w/F_s);
plot(t,real(y))
hold on
plot(t0,real(CP),'r')
set(gca, 'FontSize',14);
title('Cyclic Prefix')
xlabel('Time (sec)', 'FontSize',14);
ylabel('Amplitude', 'FontSize',14);
legend('Sine Wave', 'Cyclic Extention')
axis([0 t(end) -1 1])
