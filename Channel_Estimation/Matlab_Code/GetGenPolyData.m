function [ genPoly, n, codingRate, Memory, numStates, currentStates, statesBinary, numBranches, BranchLogic ] = GetGenPolyData
% GetGenPolyData
%
% Calculates Code Rate, Memory, States, Branches
%
% Usage :
%               [ genPoly, n, codingRate, Memory, numStates, statesDec, statesBinary, numBranches, branchesBinary ] = GetGenPolyData
%        

genPoly1 = [1 0 1 1 0 1 1];                 % Polynomials for 3GPP
genPoly2 = [1 1 1 1 0 0 1];
genPoly3 = [1 1 1 0 1 0 1];
% 
% % genPoly1 = [1 1 1];                 % Polynomials for Test
% % genPoly2 = [1 0 1];
genPoly = [genPoly1;genPoly2;genPoly3];
% genPoly = [1 0 0 0 0 1 1];

[n, polyLength] = size(genPoly);            % Number and Length of Polynomials
k = 1;                                  
codingRate = (k/n);                         % Coding Rate
Memory = polyLength-1;                      % Memory in shift register
numStates = 2^Memory;                       % Number of State for Veterbi Algorithm
statesDec = (0:numStates-1);
currentStates = statesDec+1;                % State Row Number i.e 000 is row 1
statesBinary = Dec2Bi(statesDec, Memory);   % States in Binary form
numBranches = 2^n;                          % Number of Branches per Stage in Veterbi Algorithm
[ BranchLogic ] = Dec2Bi(0:(2^n) -1, n);    % Branch Values in Binary Form
end

