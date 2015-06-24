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

% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.1
% Date		: 25 May 2008
% Details	: Corrected the 
% Version	: 1.0
% Date		: 30 April 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script for plotting the eye diagram where transmit filtering
% is performed by raised cosine filtering with alpha=0.5, alpha=1.

clear all
N  = 10^3; % number of symbols
am = 2*(rand(1,N)>0.5)-1 + j*(2*(rand(1,N)>0.5)-1); % generating random binary sequence
fs = 10; % sampling frequency in Hz

% defining the sinc filter
sincNum = sin(pi*[-fs:1/fs:fs]); % numerator of the sinc function
sincDen = (pi*[-fs:1/fs:fs]); % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

% raised cosine filter
alpha = 0.5;
cosNum = cos(alpha*pi*[-fs:1/fs:fs]);
cosDen = (1-(2*alpha*[-fs:1/fs:fs]).^2);
cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;

gt_alpha5 = sincOp.*cosOp;

alpha = 1;
cosNum = cos(alpha*pi*[-fs:1/fs:fs]);
cosDen = (1-(2*alpha*[-fs:1/fs:fs]).^2);
cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;
gt_alpha1 = sincOp.*cosOp;

% upsampling the transmit sequence 
amUpSampled = [am;zeros(fs-1,length(am))];
amU = amUpSampled(:).';

% filtered sequence
st_alpha5 = conv(amU,gt_alpha5);
st_alpha1 = conv(amU,gt_alpha1);

% taking only the first 10000 samples 
st_alpha5 = st_alpha5([1:10000]);
st_alpha1 = st_alpha1([1:10000]);

st_alpha5_reshape = reshape(st_alpha5,fs*2,N*fs/20).';
st_alpha1_reshape = reshape(st_alpha1,fs*2,N*fs/20).';

close all
figure;
plot([0:1/fs:1.99],real(st_alpha5_reshape).','b');   
title('eye diagram with alpha=0.5');
xlabel('time')
ylabel('amplitude') 
axis([0 2 -1.5 1.5])
grid on

figure;
plot([0:1/fs:1.99],real(st_alpha1_reshape).','b');  
title('eye diagram with alpha=1') 
xlabel('time')
ylabel('amplitude') 
axis([0 2 -1.5 1.5 ])
grid on