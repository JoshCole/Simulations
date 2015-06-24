function [ x, xStart, xEnd] = GetRef( RxSymbol, Nref, index )
NtxRef = round(Nref/2);
Nx = round(length(RxSymbol)/2);
xStart = Nx-NtxRef+index;
xEnd = Nx+NtxRef-1+index;

x =RxSymbol(xStart:xEnd);

end

