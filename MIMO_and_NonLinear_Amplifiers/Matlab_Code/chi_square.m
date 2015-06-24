close all
n = 100000;
sigma=1;
x=sigma*randn(1,n);
y=x.^2;
N=floor(n/100);
A=min(y); B=max(y);
Delta=(B-A)/N;
t=A-Delta/2+[1:N]*Delta;
f=hist(y,t)/(Delta*n);
subplot(2,1,1)
bar(t,f)
title('Estimated chi-square PDF')
p=Delta*f;
CDF=cumsum(p);
subplot(2,1,2);
plot(t,CDF)
title('Estimated chi-square CDF')