function [ di ] = DeInterleaver1st( ti, Fi, Ni , TTI )
% DeInterleaver1st
%
% Reasembles data in the original encoded format
%
% Usage :
%               [ codeWord, nextState ] = EncodeBranch( register, genPoly,
%               numGenPoly, Memory)
%
% Where         ti  = Is the Recived Packet, Coded or Uncoded
%               Fi  = The Number of Frames
%               Ni  = The Number of Bits per frame
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

C1 = c1Array(TTI);
R1 = Ni;
yi = zeros(R1,C1);      % Interleaver

% Create DeInterleaving Matrix
% --------------------------
%     xi1,      xi2,        xi3         ... xic1
%     xi(c1+1), xi(c1+2),   xi(c1+3)       ,xi(2*c1)
%   etc
for i = 0:C1-1
    Index = i*R1;
    for ii = 1:R1
        yi(ii,i+1) = ti(Index+ii);
    end
end

% Use inter column permutation pattern
% ------------------------------------
yi = yi(:,InterColumnPermutation+1);

% Reasemble by column
% -------------------
di = zeros(1,Ni*Fi);
for i = 0:Fi-1
    di(1+i:Fi:end) = yi(:,i+1);
end
end

