% Generate Resource Blocks
N_RB_sc = 12;
N_SC_FDMA = 7;
rb1 = zeros(N_RB_sc,N_SC_FDMA);
rb2 = zeros(N_RB_sc,N_SC_FDMA);

Method  = 'PUCCH1a';
UL_N_symb
switch Method
    case 'PUCCH1a'
        % PUCCH Format 1a
        nBits = 1;
        b = rand(1,1)>0.5;
        [ idx ] = Bi2Dec(b, nBits)+1;
        symbolArray = [1 -1];
        d = symbolArray(idx);
    case 'PUCCH1b'
        % PUCCH Format 1b
        nBits = 2;
        b = rand(1,2)>0.5;
        [ idx ] = Bi2Dec(b, nBits)+1;
        symbolArray = [1 -1i 1i -1];
        d = symbolArray(idx);
end
N_PUCCH_seq = 12;
ns = 1                    % Slot Number
l = 1;                       % Symbol Number
    c_i = 0;
    for i = 0:7
        x1 = zeros(1,31);
        x1(1) = 1;
        n=8*N_UL_sym.ns+8*l+i;
%         x1(n+31) = 
    end