function [ yi, Xi, R1] = Interleaver1st( xi,TTI )
% Interleaver1st
%
% Mixes up the packet data in accordance to the Inter Column Permutation to
% reduce the effects of burst errors
%
% Usage :
%               [ codeWord, nextState ] = EncodeBranch( register, genPoly,
%               numGenPoly, Memory)
%
% Where         xi = Is the current Packet, Coded or Uncoded
%               TTi = Transmision Time Interval Index

% Inter Column Permutation, Reduces burst error by effectily mixing the Packet up
% -------------------------------------------------------------------------------
switch TTI
case 1
   InterColumnPermutation=0;
case 2
   InterColumnPermutation=[0,1];
case 3
   InterColumnPermutation=[0,2,1,3];
case 4
   InterColumnPermutation=[0,4,2,6,1,5,3,7];
otherwise
   error('Incorrect TTI specified, must be either 1, 2, 4, 8.');
end

c1Array= [1 2 4 8];     % Number of Columns for associated TTI Index i.e Index 4 splits Packet into 8 frames

Xi = length(xi);        % Length of current Packet Coded or Uncoded
C1 = c1Array(TTI);      % Number of Frames
R1 = Xi/C1;             % Number of Rows
yi = zeros(R1,C1);      % Interleaver


% Create Interleaving Matrix
% --------------------------
%     xi1,      xi2,        xi3         ... xic1
%     xi(c1+1), xi(c1+2),   xi(c1+3)       ,xi(2*c1)
%   etc
for i = 0:R1-1
    c1Index = i*C1;
    for ii = 1:C1
        yi(i+1,ii) = xi(c1Index+ii);
    end
end

% Use inter column permutation pattern
% ------------------------------------
yi = yi(:,InterColumnPermutation+1);


end

