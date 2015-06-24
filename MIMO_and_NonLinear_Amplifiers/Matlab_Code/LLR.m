function [ LLRData ] = LLR( BitArray, RxSymbol, k, NoiseVar )
% LLR
%
% This function calculates the Maximum Log Likliehood ratio
%
% Usage :
%
% [ LLRData ] = LLR( BitArray, RxSymbol, k, NoiseVar )
%
% Where         
%               BitArray         = Symbol array for logic values [0,1] for b0,b1 ...bk
%
%				RxSymbol         = Recieved Signal
%
%				k                = Number of Bits for M-ary
%
%				NoiseVar		 = Noise Variance for system


% Size and Array definitions
% --------------------------
RxSymbolNum = length(RxSymbol);
LLRData = zeros(k,RxSymbolNum);
BitArray_0 = BitArray(:,:,1);
BitArray_1 = BitArray(:,:,2);

parfor i = 1:RxSymbolNum            % Parrallel for loop
    LLRData(:,i) =(1/(2*(NoiseVar^2)))*...
        (-(min(abs(RxSymbol(i)-BitArray_1).^2,[],2))...
        +(min(abs(RxSymbol(i)-BitArray_0).^2,[],2)));

    % Log Likehood Ratio, expontial tends to inf or -inf for high SNR due
    % to matlab resolution
    %     LLRData(:,i) = ...
%         log(sum((1/(sqrt(2*pi*NoiseVar)))*(exp(-(1/(2*(NoiseVar^2)))*((abs(RxSymbol(i)-BitArray(:,:,2)).^2)))),2) ...
%         ./
%         sum((1/(sqrt(2*pi*NoiseVar)))*(exp(-(1/(2*(NoiseVar^2)))*((abs(RxSymbol(i)-BitArray(:,:,1)).^2)))),2));
end


end

