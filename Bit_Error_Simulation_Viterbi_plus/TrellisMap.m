function [ Value, PreviousState, Decode ] = TrellisMap( Decision, BMetric, BranchNum, BlockCode, NextState, Value, PreviousState, Decode, CurrentStates, stages, states, numInput )
% TrellisMap
%
% This function Traverse through trellis tree, updating the associated
% parameters for soft decision or hard decision.
%
% Usage :
%
% [ Value, PreviousState, Decode ] = TrellisMap( Decision, BMetric,
% BranchNum, BlockCode, NextState, Value, PreviousState, Decode, CurrentStates, stages, states, numInput )
%
% Where         Decision            = HardDecision or Soft
%
%				BMetric             = Branch Metric for Hard or Soft
%                                     Decicion
%
%				BranchNum           = Branch Numbers for each state in
%                                     terms of row position
%
%				BlockCode           = Block code is either LLR for soft
%                                     Decision or binary values for Hard Decision
%
%               NextState           = Is the next state in terms of row
%                                     position
%
%               Value               = The Value of each node in the Trellis
%
%               PreviousState       = Previous state in terms of row
%                                     position
%
%               Decode              = A Log 1 or 0 depending of which
%                                     branch was chosen. Determines all
%                                     potential decoded messages
%
%               CurrentStates       = The current State in terms of row
%                                     position
%
%               stages              = Number of Stages in Trellis
%
%               states              = Number of States in Trellis
%
%               numInput            = Number is 2 for logic [0,1]

% Calculation of Branch Metrics
% -----------------------------
switch Decision
    case 'SoftDecision'
        LLRMetric = zeros(states,numInput,stages-1);
        for i = 1:stages-1
            for ii = 1:numInput
                LLRMetric(:,ii,i)= BMetric(BranchNum(:,ii),i);      % Determines the LLR in relation to the branch values for each state
            end                                                     % Branch numers refer to a row number for branch metric

        end
    case 'HardDecision'
        hammDistance = zeros(states,numInput,stages-1 );
        for i = 1:stages-1
            currentCodeWord = repmat(BlockCode(i,:),states,1);                                      % Replicates the code value for all states
            for ii = 1:numInput
                hammDistance(:,ii,i) = sum(BMetric(BranchNum(:,ii),:)~=currentCodeWord,numInput);   % Determines the Hamming Distance in relation to the branch values for each state
            end                                                                                     % Branch numers refer to a row number for branch metric
        end

end
 accMetric = zeros(1,2);
switch Decision
    case 'SoftDecision'
        for i = 1:stages-1
            for ii = 1:states
                if Value(ii,i) ~= -100000                               % Check to see if node has been initialised
                    accMetric(1) = Value(ii,i) + LLRMetric(ii,1,i);     % Cacluate accumulated metric value for Logic [0,1]
                    accMetric(2) = Value(ii,i) + LLRMetric(ii,2,i);
                    if accMetric(1) >= Value(NextState(ii,1),i+1)       % Update node if accumulated metric is less than the value of the node in question
                        Value(NextState(ii,1),i+1) = accMetric(1);
                        PreviousState(NextState(ii,1),i+1) = CurrentStates(ii);
                        Decode(NextState(ii,1),i+1) = 0;
                    end
                    if accMetric(2) >= Value(NextState(ii,2),i+1)
                        Value(NextState(ii,2),i+1) = accMetric(2);
                        PreviousState(NextState(ii,2),i+1) = CurrentStates(ii);
                        Decode(NextState(ii,2),i+1) = 1;
                    end
                end
            end
        end
    case 'HardDecision'
        for i = 1:stages-1
            for ii = 1:states
                if Value(ii,i) ~= 100000                                % Check to see if node has been initialised
                    accMetric(1) = Value(ii,i) + hammDistance(ii,1,i);  % Cacluate accumulated metric value for Logic [0,1]
                    accMetric(2) = Value(ii,i) + hammDistance(ii,2,i);  
                    if accMetric <= Value(NextState(ii,1),i+1)          % Update node if accumulated metric is less than the value of the node in question
                        Value(NextState(ii,1),i+1) = accMetric(1);
                        PreviousState(NextState(ii,1),i+1) = CurrentStates(ii);
                        Decode(NextState(ii,1),i+1) = 0;
                    end
                    if accMetric <= Value(NextState(ii,2),i+1)
                        Value(NextState(ii,2),i+1) = accMetric(2);
                        PreviousState(NextState(ii,2),i+1) = CurrentStates(ii);
                        Decode(NextState(ii,2),i+1) = 1;
                    end
                end
            end
        end
end
end
