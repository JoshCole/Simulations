function [ rateChangeArray, testArray ] = rateMatching( xi, Xi, DeltaNij,  Eini , Eplus , Eminus, rateChangeArray, arrayRows, Fi, testArray )
arrayIndex = 1;

if DeltaNij<0
    E=Eini;                 % Inital error between current and desired puncturing ratio
    m=1;                    % Index of current bit
    tempIndex = 0;
    while(m<=Xi)
        E = E-Eminus;       % Update Error
        if E<=0             % Check if bit number m should be punctured
            testArray(m,Fi+1) = NaN;
            arrayIndex = arrayIndex -1;
            E = E+Eplus;
        end
        if tempIndex ~= arrayIndex
            rateChangeArray(arrayIndex,Fi+1) = xi(m);
            tempIndex = arrayIndex;
        end
        arrayIndex = arrayIndex +1;
        m = m+1;
    end
    
else
    E = Eini;               % Initial error between current and desired puncturing ratio
    m = 1;                  % Index of current bit
    while m <=Xi
        E = E-Eminus;       % Update error
        rateChangeArray(arrayIndex,Fi+1) = xi(m);
        testArray(arrayIndex,Fi+1) = xi(m);
        while E<=0          % Check if bit number m should be repeared
            arrayIndex = arrayIndex+1;
            rateChangeArray(arrayIndex,Fi+1) = xi(m);
            testArray(arrayIndex,Fi+1) = NaN;
            E=E+Eplus;      % Update error
        end
        arrayIndex = arrayIndex+1;
        m = m+1;            % Next bit
    end
end


