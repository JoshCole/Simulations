LengthTxPac = 1000;
m =500;
    lenUpdate = m;
    if m>1
        unqArray = randi([1 LengthTxPac],[1 lenUpdate]);
        BitToChange = unique(unqArray);
        lenUpdate=m;
        while lenUpdate
            BitToChange = unique([BitToChange,unique(unqArray)]);
            lenUpdate = m -length(BitToChange);
            unqArray = randi([1 LengthTxPac],[1 lenUpdate]);
            
        end
    else
            BitToChange = randi([1 LengthTxPac],[1 m]);
    end
    unique(BitToChange)