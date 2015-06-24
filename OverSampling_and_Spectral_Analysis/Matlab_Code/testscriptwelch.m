Fs = 1000;   t = 0:1/Fs:.296;
x = cos(2*pi*t*200)+randn(size(t));  % A cosine of 200Hz plus noise
M = length(x);
win = [];
noverlap = [];

if isempty(win),
    % Use 8 sections, determine their length
    if isempty(noverlap),
        % Use 50% overlap
        L = fix(M./4.5);
        noverlap = fix(0.5.*L);
    else
        L = fix((M+7.*noverlap)./8);
    end
    % Use a default window
    win = Hanning(L);
%     win = hamming(L);
else
    % Determine the window and its length (equal to the length of the segments)
    if ~any(size(win) <= 1) | ischar(win),
        errid = generatemsgid('invalidWindow');
        errmsg = 'The WINDOW argument must be a vector or a scalar.';
        return
    elseif length(win) > 1,
        % WIN is a vector
        L = length(win);
    elseif length(win) == 1,
        L = win;
        win = Hanning(win);
    end
    if isempty(noverlap),
        % Use 50% overlap
        noverlap = fix(0.5.*L);
    end
end

% Do some argument validation
if L > M,
    errid = generatemsgid('invalidSegmentLength');
    errmsg = 'The length of the segments cannot be greater than the length of the input signal.';
    return;
end

if noverlap >= L,
    errid = generatemsgid('invalidNoverlap');
    errmsg = 'The number of samples to overlap must be less than the length of the segments.';
    return;
end
%%
%WELCH_OPTIONS   Parse the optional inputs to the PWELCH function.
%   WELCH_OPTIONS returns a structure, OPTIONS, with following fields:
%
%   options.nfft         - number of freq. points at which the psd is estimated
%   options.Fs           - sampling freq. if any
%   options.range        - 'onesided' or 'twosided' psd

% Generate defaults
options.nfft = max(256,2^nextpow2(L));
options.Fs = []; % Work in rad/sample

%%
% Compute the number of segments
k = (M-noverlap)./(L-noverlap);

% Uncomment the following line to produce a warning each time the data
% segmentation does not produce an integer number of segements.
%if fix(k) ~= k),
%   warning(generatemsgid('MustBeInteger'),'The number of segments is not an integer, truncating data.');
%end

k = fix(k);
%%
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
% range = options.range;

LminusOverlap = L-noverlap;
xStart = 1:LminusOverlap:k*LminusOverlap;
xEnd   = xStart+L-1;
esttype = 'psd';
for i = 1:k,
    [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),win,...
        options.nfft,esttype,options.Fs);
    Sxx  = Sxx + Sxxk;
end
% switch esttype,
%     case {'ms','psd'},
%         for i = 1:k,
%             [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),win,...
%                 options.nfft,esttype,options.Fs);
%             Sxx  = Sxx + Sxxk;
%         end
% 
%     case 'cpsd',
%         for i = 1:k,
%             [Sxxk,w] =  computeperiodogram({x(xStart(i):xEnd(i)),...
%                 y(xStart(i):xEnd(i))},win,options.nfft,esttype,options.Fs);
%             Sxx  = Sxx + Sxxk;
%         end
% 
%     case 'tfe'
%         Sxy = zeros(options.nfft,1); % Initialize
%         for i = 1:k,
%             [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),...
%                 win,options.nfft,esttype,options.Fs);
%             [Syxk,w] = computeperiodogram({y(xStart(i):xEnd(i)),...
%                 x(xStart(i):xEnd(i))},win,options.nfft,esttype,options.Fs);
%             Sxx  = Sxx + Sxxk;          
%             Sxy  = Sxy + Syxk;
%         end
% 
%     case 'mscohere'
%         % Note: (Sxy1+Sxy2)/(Sxx1+Sxx2) != (Sxy1/Sxy2) + (Sxx1/Sxx2)
%         % ie, we can't push the computation of Cxy into computeperiodogram.
%         Sxy = zeros(options.nfft,1); % Initialize
%         Syy = zeros(options.nfft,1); % Initialize
%         for i = 1:k,
%             [Sxxk,w] = computeperiodogram(x(xStart(i):xEnd(i)),...
%                 win,options.nfft,esttype,options.Fs);
%             [Syyk,w] =  computeperiodogram(y(xStart(i):xEnd(i)),...
%                 win,options.nfft,esttype,options.Fs);
%             [Sxyk,w] = computeperiodogram({x(xStart(i):xEnd(i)),...
%                 y(xStart(i):xEnd(i))},win,options.nfft,esttype,options.Fs);
%             Sxx  = Sxx + Sxxk;
%             Syy  = Syy + Syyk;
%             Sxy  = Sxy + Sxyk;
%         end
% end
Sxx = Sxx./k; % Average the sum of the periodograms

% if any(strcmpi(esttype,{'tfe','mscohere'})),
%     Sxy = Sxy./k; % Average the sum of the periodograms
% 
%     if strcmpi(esttype,'mscohere'),
%         Syy = Syy./k; % Average the sum of the periodograms
%     end
% end