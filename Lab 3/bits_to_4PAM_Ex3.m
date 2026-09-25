function symbols = bits_to_4PAM_Ex3(bitseq,A)
symbols= zeros(size(bitseq(1:end/2)));

for k=0 : length(symbols)-1
    if(bitseq(2*k+1)==0 && bitseq(2*k+2)== 0)
        symbols(k+1) = 3*A;
    elseif(bitseq(2*k+1)==0 && bitseq(2*k+2)== 1)
        symbols(k+1) = 1*A;
    elseif(bitseq(2*k+1)==1 && bitseq(2*k+2)== 1)
        symbols(k+1) = -1*A;
    elseif(bitseq(2*k+1)==1 && bitseq(2*k+2)== 0)
        symbols(k+1) = -3*A;
    end 
end

end