function [symbols] = bits_to_4PAM(bits1, bits2)
symbols= zeros(size(bits1));

for k=1: length(bits1)
    if(bits1(k)==0 && bits2(k)== 0)
        symbols(k)= 3;
    elseif(bits1(k)==0 && bits2(k)== 1)
        symbols(k)= 1;
    elseif(bits1(k)==1 && bits2(k)== 1)
        symbols(k)= -1;
    elseif(bits1(k)==1 && bits2(k)== 0)
        symbols(k)= -3;
    end 
end

end