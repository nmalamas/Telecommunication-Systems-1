function [X] = bits_to_2PAM(bit_seq)
X = zeros(1, length(bit_seq));

for k=1:length(bit_seq)
    if(bit_seq(k) == 0)
        X(k) = 1;
    elseif(bit_seq(k) == 1)
        X(k) = -1;
    end
end
end