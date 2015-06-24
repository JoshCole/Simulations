function [P,f] = computeperiodogram(x,win,nfft,esttype,varargin)
%COMPUTEPERIODOGRAM   Periodogram spectral estimation.
%   This function is used to calculate the Power Spectrum Sxx, and the
%   Cross Power Spectrum Sxy.
%
%   Sxx = COMPUTEPERIODOGRAM(X,WIN,NFFT) where x is a vector returns the
%   Power Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Sxy = COMPUTEPERIODOGRAM({X,Y},WIN,NFFT) returns the Cross Power
%   Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Inputs:
%    X           - Signal vector or a cell array of two elements containing
%                  two signal vectors.
%    WIN         - Window
%    NFFT        - Number of frequency points (FFT) or vector of
%    frequencies at whicch periodogram is desired (Goertzel)
%    WINCOMPFLAG - A string indicating the type of window compensation to
%                  be done. The choices are: 
%                  'ms'  - compensate for Mean-square (Power) Spectrum;
%                          maintain the correct power peak heights.
%                  'psd' - compensate for Power Spectral Denstiy (PSD);
%                          maintain correct area under the PSD curve.
%
%   Output:
%    Sxx         - Power spectrum [Power] over the whole Nyquist interval. 
%      or
%    Sxy         - Cross power spectrum [Power] over the whole Nyquist
%                  interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2007/12/14 15:15:02 $ 

error(nargchk(3,5,nargin,'struct'));
if nargin < 4,
    esttype = 'psd'; % Default, compenstate for window's power.
end

if nargin < 5 || isempty(varargin{1}),
    Fs = 2*pi;
else
    Fs = varargin{1};     
end

% Validate inputs and convert row vectors to column vectors.
[x,Lx,y,is2sig,win,errid,errmsg] = validateinputs(x,win,nfft);
if ~isempty(errmsg), error(errid, errmsg); end

% Determine if FFT or Goertzel
freqVectorSpecified = false;
if (length(nfft)> 1), freqVectorSpecified = true; end

% Window the data
xw = x.*win;
if is2sig, yw = y.*win; end 

% Evaluate the window normalization constant.  A 1/N factor has been
% omitted since it will cancel below.
if strcmpi(esttype,'ms'),
    % The window is convolved with every power spectrum peak, therefore
    % compensate for the DC value squared to obtain correct peak heights.
    U = sum(win)^2;
else
    U = win'*win;  % compensates for the power of the window.
end

% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels

[Xx,f] = computeDFT(xw,nfft,Fs);
if is2sig, [Yy,f] = computeDFT(yw,nfft,Fs); end

P = Xx.*conj(Xx)/U;      % Auto spectrum.
if is2sig,
    P = Xx.*conj(Yy)/U;  % Cross spectrum.
end

%--------------------------------------------------------------------------
function [x,Lx,y,is2sig,win,errid,errmsg] = validateinputs(x,win,nfft)
% Validate the inputs to computexperiodogram and convert the input signal
% vectors into column vectors.

errid = '';
errmsg = ''; % in case of early return.

% Set defaults and convert to row vectors to columns.
y     = [];
Ly    = 0;
is2sig= false;
win   = win(:);
Lw    = length(win);

% Determine if one or two signal vectors was specified.
Lx = length(x);
if iscell(x),
    if length(x) > 1,
        y = x{2};
        y = y(:);
        is2sig = true;
    end
    x = x{1};
    Lx = length(x);
end
x = x(:);

if is2sig,
    Ly  = length(y);
    if Lx ~= Ly,
        errid = generatemsgid('invalidInputSignalLength');
        errmsg = 'The length of the two input vectors must be equal to calculate the cross spectral density.';
        return;
    end
end

if Lx ~= Lw,
    errid = generatemsgid('invalidWindow');
    errmsg = 'The WINDOW must be a vector of the same length as the signal.';
    return;
end


% [EOF]
