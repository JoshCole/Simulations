function [bb, tim] = rcosfir(r, N_T, rate, T, fil_type, col) 
%RCOSFIR Designs raised cosine FIR filter. 
%       NOTE: The output of this function is a raised cosine FIR filter. You 
%             should use function RCOSFLT to filter your digital signal 
%             instead of using FILTER. 
%       RCOSFIR(R, N_T, RATE, T) produces time response and frequency response 
%       of a raised cosine filter at given rolloff factor R. N_T is a length 2 
%       vector that indicates the number of T windows in the design. The 
%       default value for N_T is [-3, 3]. RATE is the sample point in each T, 
%       or, the sampling rate of the filter is T/RATE. The default value of  
%       RATE is 5. The order of the FIR filter is determined by  
%           (N_T(2) - N_T(1)) * RATE 
%       T is the symbol interval. The default value of T is 1. 
%       The time response of the raised cosine filter has the form of 
%       h(t) = sinc(t/T) cos(pi R t/T)/(1 - 4 R^2 t^2 /T^2) 
%       The frequency domain has the spectrum as  
%              /  T                                 when 0 < |f| < (1-r)/2/T 
%              |          pi T         1-R    T           1-R         1+R 
%       H(f) = < (1 + cos(----) (|f| - ----) ---    when  --- < |f| < --- 
%              |            r           2T    2           2 T         2 T 
%              \  0                                 when |f| > (1+r)/2/T 
% 
% 
%       RCOSFIR(R, N_T, RATE, T, FILTER_TYPE) produces time response and 
%       frequency response of a square root raised cosine filter if  
%       FILTER_TYPE == 'sqrt'. 
% 
%       RCOSFIR(R, N_T, RATE, T, FILTER_TYPE, COLOR) produces time response 
%       and frequency response with the curve color as specified in the string 
%       variable COLOR. The string in COLOR can be any type as defined in 
%       PLOT. 
% 
%       B = RCOSFIR(...) returns the designed raised cosine FIR filter. 
%       The size of vector B is  (N_T(2) - N_T(1))*RATE + 1.  
% 
%       [B, Sample_Time] = RCOSFIR(...) returns the FIR filter and the sample 
%       time for the filter. Note that the filter sample time is T / RATE. 
%       The digital transmission rate is T. 
% 
%       See also RCOSIIR, RCOSFLT. 
 
%       Wes Wang 5/25/94, 10/11/95 
%       Copyright (c) 1995-96 by The MathWorks, Inc. 
%       $Revision: 1.1 $  $Date: 1996/04/01 18:02:13 $ 
 
%routine check 
if nargin < 1 
    error('Not enough input variable for RCOSFIR') 
elseif nargin < 2 
    N_T = [3 3]; rate = 5; T = 1; fil_type = 'normal'; 
elseif nargin < 3,  
    rate = 5; T = 1; fil_type = 'normal'; 
elseif nargin < 4,  
    T = 1; fil_type = 'normal'; 
elseif nargin < 5,  
    fil_type = 'normal'; 
end; 
 
if r < 0,  
    r = 0; 
elseif r > 1, 
    r = 1; 
end; 
 
[N_T, rate, T, fil_type] = checkinp(N_T, rate, T, fil_type,... 
                                    [3 3], 5,  1, 'normal'); 
if length(N_T) < 2 
    N_T = [N_T N_T]; 
end; 
 
if (rate < 1) | (ceil(rate) ~= rate) 
    error('RATE in RCOSFIR must be a positive integer') 
end 
 
% calculation 
N_T(1) = -abs(N_T(1)); 
%time_T = [N_T(1) : 1/rate : N_T(2)];  
time_T = [0 : 1/rate : max(N_T(2), abs(N_T(1)))];  
cal_time = time_T * T; 
time_T_r = r * time_T; 
if ~isempty(findstr(fil_type,'root')) | ... 
        ~isempty(findstr(fil_type,'sqrt'))%square root 
    if r==0  %use rcosdev we have h(t) = sinc(pi*t/T) / sqrt(T) 
        b = sinc(time_T) /sqrt(T); 
    else 
        den = 1 - (4 * time_T_r).^2; 
        ind1 = find(den ~= 0); 
        ind2 = find(den == 0); 
 
        % at the regular point use the regular rule 
        nterm1 = cos((1 + r) * pi * time_T); 
        % term 2 in the numerator: 
        nterm2 = sinc((1-r)*time_T) *((1-r) * pi / 4 / r); 
        b(ind1) = (nterm1(ind1) + nterm2(ind1)) ./ den(ind1) * 4 * r / pi / sqrt(T); 
 
        % Note that time_T_r may not be zero at the point of den is zero. 
        % Use L'Hopital rule to compute the value at those points. 
        if ~isempty(ind2) 
            nterm1 = (1 + r) * pi * sin((1 + r) * pi * time_T(ind2)); 
            nterm2 = (sin((1-r) * pi * time_T(ind2)) -... 
                 (1 - r) * pi * time_T(ind2) .* cos((1 - r) * pi * time_T(ind2))) / ... 
                  4 / r ./ time_T(ind2) .^ 2; 
            b(ind2) = (nterm1 + nterm2) ./ time_T_r(ind2) / 8 / sqrt(T) / pi; 
        end; 
    end; 
%    b = b * sqrt(T); 
    b = b * sqrt(cal_time(2)); 
else    % regular case 
    den = 1 - (2 * time_T_r).^2; 
    ind1 = find(den ~= 0); 
    ind2 = find(den == 0); 
 
    % At the regular points, use the regular rule 
    b(ind1) = sinc(time_T(ind1)) .* ... 
        cos(pi * time_T_r(ind1)) ./ ... 
            den(ind1); 
 
    % At the points denominator equals to zero, use L'Hopital rule 
    % Note that time_T_r will not be zero at the point of den is zero. 
    if ~isempty(ind2) 
         b(ind2) = pi * sinc(time_T(ind2)) .* ... 
             sin(pi * time_T_r(ind2)) ./ ... 
                 time_T_r(ind2) / 8; 
    end; 
end; 
 
b = [b(rate * abs(N_T(1))+1 : -1 : 1), b(2 : rate * N_T(2)+1)]; 
tim = cal_time(2) - cal_time(1); 
 
% In the case needs a plot 
if nargout < 1 
    if nargin < 6 
        col = 'y-'; 
    end; 
 
    % the time response part 
    hand = subplot(211); 
    out = filter(b, 1, [1, zeros(1, length(cal_time) - 1)]); 
    plot(cal_time, out, col) 
    % if not hold, change the axes 
    hol = get(hand,'NextPlot'); 
    if (hol(1:2) ~= 'ad') | (max(get(hand,'Ylim')) < max(b)) 
        axis([min(cal_time), max(cal_time), min(out) * 1.1, max(out) * 1.1]); 
        xlabel('time'); 
        title('Impulse Response of the Raised Cosine Filter (with time shift)') 
    end; 
 
    % the frequency response part 
    hand = subplot(212); 
    len = length(b); 
    P = abs(fft(b)) * abs(N_T(2) - N_T(1)) / len * T; 
    f = (0 : len / 2) / len * rate / T; 
    ind = find(f < 1.5 / T); 
    f = f(ind); 
    P = P(ind); 
    plot(f, P, col); 
    hol = get(hand, 'NextPlot'); 
    if hol(1:2) ~= 'ad' 
        xlabel('frequency'); 
        ylabel('Amplitude'); 
        title('Frequency Response of the Raised Cosine Filter') 
    end; 
else  
    bb = b; 
end; 
%--end of rcosfir.m-- 