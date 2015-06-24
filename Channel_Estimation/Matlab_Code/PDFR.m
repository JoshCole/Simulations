close all
clear all
N = 10^6;
x = randn(1,N); % gaussian random variable, mean 0, variance 1
y = randn(1,N); % gaussian random variable, mean 0, variance 1
z = (x + j*y);  % complex random variable

% probability density function of abs(z)
zBin = [0:0.01:7];
sigma2 = 1;
pzTheory = (zBin/sigma2).*exp(-(zBin.^2)/(2*sigma2)); % theory
[nzSim zBinSim] = hist(abs(z),zBin); % simulation

% probability density of theta
thetaBin = [-pi:0.01:pi];
pThetaTheory = 1/(2*pi)*ones(size(thetaBin));
[nThetaSim thetaBinSim] = hist(angle(z),thetaBin); % simulation

figure
plot(zBinSim,nzSim/(N*0.01),'m','LineWidth',2);
hold on
plot(zBin,pzTheory,'b.-')
xlabel('z');
ylabel('probability density, p(z)');
legend('simulation','theory');
title('Probability density function of abs(z)' )
axis([0 7 0 0.7]);
grid on
temp = nThetaSim/(N*0.01);
figure
plot(thetaBinSim,temp,'m','LineWidth',2);
hold on
plot(thetaBin,pThetaTheory,'b.-')
xlabel('theta');
ylabel('probability density, p(theta)');
legend('simulation','theory');
title('Probability density function of theta')
axis([-pi pi 0 0.2])
grid on