function varargout = welch(x,esttype,varargin)
%WELCH Welch spectral estimation method.
%   [Pxx,F] = WELCH(X,WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,ESTTYPE)
%   [Pxx,F] = WELCH({X},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'psd')
%   [Pxx,F] = WELCH({X},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'ms')
%   [Pxy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'cpsd')
%   [Txy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'tfe')
%   [Cxy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'mscohere')
%
%   Inputs:
%      see "help pwelch" for complete description of all input arguments.
%      ESTTYPE - is a string specifying the type of estimate to return, the
%                choices are: psd, cpsd, tfe, and mscohere.
%
%   Outputs:
%      Depends on the input string ESTTYPE:
%      Pxx - Power Spectral Density (PSD) estimate, or
%      MS  - Mean-square spectrum, or
%      Pxy - Cross Power Spectral Density (CPSD) estimate, or
%      Txy - Transfer Function Estimate (TFE), or
%      Cxy - Magnitude Squared Coherence.
%      F   - frequency vector, in Hz if Fs is specified, otherwise it has
%            units of rad/sample

%   Author(s): P. Pacheco
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2007/12/14 15:15:21 $


%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and
%         Modeling, John Wiley & Sons, 1996.

error(nargchk(2,8,nargin,'struct'));
error(nargoutchk(0,3,nargout,'struct'));

% Parse input arguments.
[x,M,isreal_x,y,Ly,win,winName,winParam,noverlap,k,L,options] = ...
    welchparse(x,esttype,varargin{:});

% Frequency vector was specified, return and plot two-sided PSD
freqVectorSpecified = false; nrow = 1;
if length(options.nfft) > 1, 
    freqVectorSpecified = true; 
    [ncol,nrow] = size(options.nfft); 
end

% Compute the periodogram power spectrum of each segment and average always
% compute the whole power spectrum, we force Fs = 1 to get a PS not a PSD.

% Initialize
if freqVectorSpecified, 
    Sxx = zeros(length(options.nfft),1);   
else
    Sxx = zeros(options.nfft,1); 
end
range = options.range;

LminusOverlap = L-noverlap;
xStart = 1:LminusOverlap:k*LminusOverlap;
xEnd   = xStart+L-1;
switch esttype,
    case {'ms','psd'},
        for i = 1:k,
            [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),win,...
                options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
        end

    case 'cpsd',
        for i = 1:k,
            [Sxxk,w] =  computeperiodogram({x(xStart(i):xEnd(i)),...
                y(xStart(i):xEnd(i))},win,options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
        end

    case 'tfe'
        Sxy = zeros(options.nfft,1); % Initialize
        for i = 1:k,
            [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),...
                win,options.nfft,esttype,options.Fs);
            [Syxk,w] = computeperiodogram({y(xStart(i):xEnd(i)),...
                x(xStart(i):xEnd(i))},win,options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;          
            Sxy  = Sxy + Syxk;
        end

    case 'mscohere'
        % Note: (Sxy1+Sxy2)/(Sxx1+Sxx2) != (Sxy1/Sxy2) + (Sxx1/Sxx2)
        % ie, we can't push the computation of Cxy into computeperiodogram.
        Sxy = zeros(options.nfft,1); % Initialize
        Syy = zeros(options.nfft,1); % Initialize
        for i = 1:k,
            [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),...
                win,options.nfft,esttype,options.Fs);
            [Syyk,w] =  computeperiodogram(y(xStart(i):xEnd(i)),...
                win,options.nfft,esttype,options.Fs);
            [Sxyk,w] = computeperiodogram({x(xStart(i):xEnd(i)),...
                y(xStart(i):xEnd(i))},win,options.nfft,esttype,options.Fs);
            Sxx  = Sxx + Sxxk;
            Syy  = Syy + Syyk;
            Sxy  = Sxy + Sxyk;
        end
end
Sxx = Sxx./k; % Average the sum of the periodograms

if any(strcmpi(esttype,{'tfe','mscohere'})),
    Sxy = Sxy./k; % Average the sum of the periodograms

    if strcmpi(esttype,'mscohere'),
        Syy = Syy./k; % Average the sum of the periodograms
    end
end

% Generate the freq vector directly in Hz to avoid roundoff errors due to
% conversions later.
if ~freqVectorSpecified, 
    w = psdfreqvec('npts',options.nfft, 'Fs',options.Fs);
else
    if strcmpi(options.range,'onesided')
        warning(generatemsgid('InconsistentRangeOption'),...
            'Ignoring ''onesided'' option. When a frequency vector is specified, a ''twosided'' PSD is computed');
    end
    options.range = 'twosided';
end


% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, corresponding freq vector and freq units.
[Pxx,w,units] = computepsd(Sxx,w,options.range,options.nfft,options.Fs,esttype);

if any(strcmpi(esttype,{'tfe','mscohere'})),
    % Cross PSD.  The frequency vector and xunits are not used.
    [Pxy,f,xunits] = computepsd(Sxy,w,options.range,options.nfft,options.Fs,esttype);

    % Transfer function estimate.
    if strcmpi(esttype,'tfe'),
        Pxx = Pxy ./ Pxx; % Txy
    end

    % Magnitude Square Coherence estimate.
    if strcmpi(esttype,'mscohere'),
        % Auto PSD for 2nd input vector. The freq vector & xunits are not
        % used.
        [Pyy,f,xunits] = computepsd(Syy,w,options.range,options.nfft,options.Fs,esttype);
        Pxx = (abs(Pxy).^2)./(Pxx.*Pyy); % Cxy
    end
end

if nargout==0
    w = {w};
    if strcmpi(units,'Hz'), w = {w{:},'Fs',options.Fs};  end
    % Create a spectrum object to store in the Data object's metadata.
    percOverlap = (noverlap/L)*100;
    hspec = spectrum.welch({winName,winParam},L,percOverlap);

    switch lower(esttype)
        case 'tfe'
            if strcmpi(options.range,'onesided'), range='half'; else range='whole'; end
            h = dspdata.freqz(Pxx,w{:},'SpectrumRange',range);
        case 'mscohere'
            if strcmpi(options.range,'onesided'), range='half'; else range='whole'; end
            h = dspdata.magnitude(Pxx,w{:},'SpectrumRange',range);
        case 'cpsd'
            h = dspdata.cpsd(Pxx,w{:},'SpectrumType',options.range);
        otherwise
            h = dspdata.psd(Pxx,w{:},'SpectrumType',options.range);
    end
    h.Metadata.setsourcespectrum(hspec);
    plot(h);
else
   % If the frequency vector was specified as a row vector, return outputs 
   % the correct dimensions
   if nrow > 1,  
       Pxx = Pxx.'; w = w.'; 
   end
    varargout = {Pxx,w}; % Pxx=PSD, MEANSQUARE, CPSD, or TFE
end

% [EOF]
