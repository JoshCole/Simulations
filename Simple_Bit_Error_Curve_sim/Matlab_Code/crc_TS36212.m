function [ A, B ] = crc_TS36212(mode, packet, parity )
% Parity configuration
% --------------------
switch parity
    case 'gCRC24A'
        gCRC = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1];
    case 'gCRC24B'
        gCRC = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1];
    case 'gCRC16'
        gCRC = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
    case 'gCRC8'
        gCRC = [1 1 0 0 1 1 0 1 1];
end

% Frame adjustment in respect to Tx or Rx
% ---------------------------------------
if strcmp(mode,'Transmitter') == true
    [M N]=size(gCRC);
    newpacket = [packet false(1,N-1)];
elseif strcmp(mode,'Reciever') == true
    newpacket = reshape(packet,1,[]);
end

% Deconvolution and polynomial division
% -------------------------------------
[q rem]=deconv(newpacket,gCRC);
rem=abs(rem);
for i=1:length(rem)
    if ( mod(rem(i),2)== 0 )
        rem(i)=false;
    else
        rem(i)=true;
    end
end

% Frame to Transmit
% -----------------
frame = bitor(newpacket,rem);
[M,newFrameSize] = size(frame);
if strcmp(mode,'Transmitter') == true
    A = frame;
    B = newFrameSize;
    
% CRC Error check
% -----------------
elseif strcmp(mode,'Reciever') == true
    frameError = false;
    A = frameError;
    B = newFrameSize;
    if sum(rem)>0
        frameError = true;
        A = frameError;
    end
end
end

