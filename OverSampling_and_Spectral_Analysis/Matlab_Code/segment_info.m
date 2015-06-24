function [L,noverlap,win,msg] = segment_info(M,win,noverlap)
%SEGMENT_INFO   Determine the information necessary to segment the input data.         
%
%   Inputs:
%      M        - An integer containing the length of the data to be segmented
%      WIN      - A scalar or vector containing the length of the window or the window respectively                 
%                 (Note that the length of the window determines the length of the segments)
%      NOVERLAP - An integer containing the number of samples to overlap (may be empty)
%
%   Outputs:
%      L        - An integer containing the length of the segments
%      NOVERLAP - An integer containing the number of samples to overlap
%      WIN      - A vector containing the window to be applied to each section
%      MSG      - A string containing possible error messages
%
%
%   The key to this function is the following equation:
%
%      K = (M-NOVERLAP)/(L-NOVERLAP)
%
%   where
%
%      K        - Number of segments
%      M        - Length of the input data X
%      NOVERLAP - Desired overlap
%      L        - Length of the segments
%   
%   The segmentation of X is based on the fact that we always know M and two of the set
%   {K,NOVERLAP,L}, hence determining the unknown quantity is trivial from the above
%   formula.

% Initialize outputs
L = [];
msg = '';

% Check that noverlap is a scalar
if any(size(noverlap) > 1),
   msg = 'You must specify an integer number of samples to overlap.';
   return
end

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
   win = hamming(L);
else
   % Determine the window and its length (equal to the length of the segments)
   if ~any(size(win) <= 1) | ischar(win),
      msg = 'The WINDOW argument must be a vector or a scalar.';
      return
   elseif length(win) > 1,
      % WIN is a vector
      L = length(win);
   elseif length(win) == 1,
      L = win;
      win = hamming(win);
   end
   if isempty(noverlap),
      % Use 50% overlap
      noverlap = fix(0.5.*L);
   end
end

% Do some argument validation
if L > M,
   msg = 'The length of the segments cannot be greater than the length of the input signal.';
   return
end

if noverlap >= L,
   msg = 'The number of samples to overlap must be less than the length of the segments.';
   return
end