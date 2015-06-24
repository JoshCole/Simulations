%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creative Commons
% Attribution-Noncommercial 2.5 India
% You are free:
% to Share — to copy, distribute and transmit the work
% to Remix — to adapt the work
% Under the following conditions:
% Attribution. You must attribute the work in the manner 
% specified by the author or licensor (but not in any way 
% that suggests that they endorse you or your use of the work). 
% Noncommercial. You may not use this work for commercial purposes. 
% For any reuse or distribution, you must make clear to others the 
% license terms of this work. The best way to do this is with a 
% link to this web page.
% Any of the above conditions can be waived if you get permission 
% from the copyright holder.
% Nothing in this license impairs or restricts the author's moral rights.
% http://creativecommons.org/licenses/by-nc/2.5/in/

% Script for simulating QPSK(4-QAM) transmission and reception and compare the 
% simulated and theoretical symbol error probability

% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 16 February 2007
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% symbol error rate for QPSK(4-QAM) modulation

clear
N = 10^5; % number of symbols
Es_N0_dB = [-3:20]; % multiple Eb/N0 values
ipHat = zeros(1,N);
for ii = 1:length(Es_N0_dB)
ip = (2*(rand(1,N)>0.5)-1) + j*(2*(rand(1,N)>0.5)-1); %
s = (1/sqrt(2))*ip; % normalization of energy to 1
n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

% demodulation
y_re = real(y); % real
y_im = imag(y); % imaginary
ipHat(find(y_re < 0 & y_im < 0)) = -1 + -1*j;
ipHat(find(y_re >= 0 & y_im > 0)) = 1 + 1*j;
ipHat(find(y_re < 0 & y_im >= 0)) = -1 + 1*j;
ipHat(find(y_re >= 0 & y_im < 0)) = 1 - 1*j;

nErr(ii) = size(find([ip- ipHat]),2); % couting the number of errors
end

simSer_QPSK = nErr/N;
theorySer_QPSK = erfc(sqrt(0.5*(10.^(Es_N0_dB/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(Es_N0_dB/10))))).^2;

close all
figure
semilogy(Es_N0_dB,theorySer_QPSK,'b.-');
hold on
semilogy(Es_N0_dB,simSer_QPSK,'mx-');
axis([-3 15 10^-5 1])
grid on
legend('theory-QPSK', 'simulation-QPSK');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for QPSK(4-QAM)')