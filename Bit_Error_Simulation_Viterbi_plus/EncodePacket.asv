function [ codeWord ] = EncodePacket( Packet, genPoly, numGenPoly, Memory, Method )
% EncodePacket
%
% This function Encodes the packet in accordance to the coding rate
%
% Usage :
%
% [ codeWord ] = EncodePacket( Packet, genPoly, numGenPoly, Method )
%
% Where         Packet      = Packet to be encoded
%
%				genPoly     = Polynomials used for encoding
%
%				numGenPoly  = Number of Polynomials (k/n)
%
%				Method      = Uncoded doesn't apply coding, Coded applies
%                             coding


if strcmp(Method,'Uncoded')     % Determines whether coding is to be applied
    codeWord = Packet;
    return
elseif strcmp(Method,'Coded')
    Packet = [Packet,zeros(1,Memory)];  % Tail Biting
    N = length(Packet);
    code = zeros(numGenPoly,N);
    
    
    for i = 1:numGenPoly
        temp = mod(conv(Packet,genPoly(i,:)),2);    % Uses Conv to determine code word
        code(i,:) = temp(1:N);                      % Extracts relevent data for Conv
    end
    codeWord = zeros(1,numGenPoly*N);
    % Assemble code word from n outputs
    % ---------------------------------
    for i = 0:numGenPoly-1
        codeWord(1+i:numGenPoly:end) = code(i+1,:);
    end
end
end

