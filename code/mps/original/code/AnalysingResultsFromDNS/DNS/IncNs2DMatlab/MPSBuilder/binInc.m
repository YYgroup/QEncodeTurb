function incBinNum  = binInc(binNum)
%Increments any given binary number by 1. The inputted binary number must
%be in the form of a vector, and must only consist of 1s and 0s. Overflow
%is handled by setting everything to zero -- length of binNum is not
%incremented.

L = length(binNum);
incComplete = 0;
i = 1;

while(incComplete == 0 && i<=L)
    %If trying to increment at a location with 1, set to 0 and move to
    %next location (i.e. "carry" it with you). When first 0 is found,
    %set it to 1. 
    if( binNum(i) == 0)
        binNum(i) = 1;
        incComplete = 1;
    else
        binNum(i) = 0;
    end
    i = i+1; 
end

incBinNum = binNum;

end