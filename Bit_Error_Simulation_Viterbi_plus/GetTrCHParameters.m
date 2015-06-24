function [ yi, DeltaNij , Nij , Fi, TrCHArray ] = GetTrCHParameters( yi, txPacketLength )
% GetTrCHParameters
%
% Determines:   Frame size and number of bits per frame
%               Whether Puncturing or repatition is required
%               Creates a new array of required size for transmission
% Usage :
%               [ DeltaNij , Nij , Fi, TrCHArray ] = GetTrCHParameters( yi,
%               txCRCPacketLength )
%
% Where         yi = Packet after 1st Interleaver
%               txCRCPacketLength = Transmitted Packet Length

[Nij,Fi] = size(yi);     % [Number of radio frames, number of bits per frame]

% Determines if Puncturing or repatition is required
% --------------------------------------------------
newRow = txPacketLength/Fi;
if (newRow)<Nij
    DeltaNij = newRow - Nij;
elseif (newRow)>Nij 
    DeltaNij = newRow -Nij;
else
    DeltaNij = 0;
end

% Calcuates transmistion array size
% ---------------------------------
Ni = Nij+DeltaNij;
TrCHArray = zeros(Ni,Fi);
yi = reshape(yi,1,[]);
end

