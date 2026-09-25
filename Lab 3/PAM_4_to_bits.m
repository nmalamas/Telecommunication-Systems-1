function est_bit = PAM_4_to_bits(X, A)
est_bit = zeros(1,length(X));

for k=0:length(X)-1
    if X(k+1) == 3*A
        est_bit(2*k+1) = 0;
        est_bit(2*k+2) = 0;
    elseif X(k+1) == 1*A
        est_bit(2*k+1) = 0;
        est_bit(2*k+2) = 1;
    elseif X(k+1) == -1*A
        est_bit(2*k+1) = 1;
        est_bit(2*k+2) = 1;
    elseif X(k+1) == -3*A
        est_bit(2*k+1) = 1;
        est_bit(2*k+2) = 0;
    end
end
