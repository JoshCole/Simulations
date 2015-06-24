function [ SoftDecision, HardDecision ] = DecisionType( LLRData, k, TrCHFrameSize )
% DecisionType
%
% This function Provides the values for Hard and Soft Decision
%
% Usage :
%
% [ SoftDecision, HardDecision ] = DecisionType( LLRData, k, TrCHFrameSize )
%
% Where         LLRData         = LLR for b0,b1 ...bk
%
%				k               = Number of Bits for M-ary
%
%				TrCHFrameSize   = Noise Variance for system
% LLRData = LLRData';
HardDecisionArray = zeros(k,TrCHFrameSize);
for i = 1:TrCHFrameSize
    for ii =1:k
        if LLRData(ii,i)>=0
            HardDecisionArray(ii,i)=1;
        end
    end
end
% Reshapes data into datastream
% -----------------------------
HardDecision = reshape(HardDecisionArray,1,[]);      
SoftDecision = reshape(LLRData,1,[]);

end

