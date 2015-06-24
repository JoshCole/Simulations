function p=pulsegen(fs,T,edge,type,f,opt)
%p=pulsegen(fs,T,edge,type,f,opt);
%a signal generation program
%fs is the sampling frequency
%T is the total signal length
%edge is a decay parameter for some waveforms
%   it is used in 'gaussian', 'monocycle',  'biexponential', 'mexican hat', 'sinc', 'double sinc', 'sinc squared'
%   and windowed sweep
%   it is mostly a parameter to describe how much the edge of the pulse is decayed.
%type is the type of the waveform desired
%   allowable types are 'gaussian', 'square', 'triangle', 'monocycle',
%   'biexponential', 'mexican hat', 'sinc', 'double sinc', 'sinc squared','sweep', and 'windowed sweep'
%f is the modulation frequency if left out it is assumed 0.
%opt is an optional argument for pulse waveforms requiring a lower and higher frequency
%    it is used in 'double sinc' ,'sweep' and 'windowed sweep' for the low and high frequency.
% the pulses are always normalized to a peak amplitude of 1

if nargin<4 %test for optional arguments
    error('not enough input arguments');
elseif nargin==4
    f=0;
    opt=[16*edge/(5*T),64*edge/(5*T)];
elseif nargin==5
    opt=[16*edge/(5*T),64*edge/(5*T)];
end
if (edge==0)
    edge==1;
end
t=-T/2:1/fs:T/2;
sig=(T/8/edge)^2;
switch type
case {'guassian'}  %generate a guassian pulse
    y=exp(-(t).^2/sig);
case {'square'}  %generate a square pulse
    y=ones(size(t));
case {'triangle'} %generate a triangle pulse
    y=(t+T/2).*(t<0)-(t-T/2).*(t>=0);
case {'monocycle'} %generate a gaussian monocycle
    y=2*t./sig.*exp(-(t).^2/sig); 
case {'biexponential'} %generate a biexponential pulse
    y=exp(-abs(t)*8*edge/T);
case {'mexican hat'}  %generate a gaussian second deriviative
    z=t./sqrt(0.75*sig);
    y=sqrt(1/2*pi).*(1-z.^2).*exp(-z.^2/2);
case {'sinc'} %generate a sinc function
    y=sinc(2*pi*edge*16.*t/(5*T));
case {'double sinc'} %generate a bandlimited function from two sinc functions
    y=opt(1)*sinc(2*opt(1).*t)-opt(2)*sinc(2*opt(2).*t);
case {'sinc squared'} %generates sinc squared function
    y=sinc(2*pi*edge*16.*t/(5*T)).^2;
case {'sweep'} %generate frequency sweep
    theta=(opt(1)+(opt(2)-opt(1))/T).*(t+T/2);
    y=real(exp(j*(theta.*(t+T/2)-pi/2)));
case {'windowed sweep'} %generate a windows frequency sweep
    theta=(opt(1)+(opt(2)-opt(1))/T).*(t+T/2);
    y=real(exp(j*(theta.*(t+T/2)-pi/2)));
    c=length(y);
    edge=min(1,edge);
    edge=max(0,edge);
    w=hamming(ceil(c*(1-edge)));
    w2=[w(1:ceil(length(w)/2));ones(c-length(w),1);w(ceil(length(w)/2)+1:end)]';
    y=w2.*y;
otherwise
    error('invalid pulse type');
end

%apply a modulation
if f~=0
    y=y.*cos(2*pi*t*f);
end
%normalize the peak of the pulse to 1
p=y./max(abs(y));
