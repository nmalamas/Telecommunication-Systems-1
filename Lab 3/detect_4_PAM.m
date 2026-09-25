function est_X = detect_4_PAM(Y, A)

est_X = zeros(1,length(Y));

for k=1:length(Y)
    if (Y(k) <= -2*A)
        est_X(k) = -3*A;
    elseif (-2*A < Y(k) && Y(k) <= 0)
        est_X(k) = -1*A;
    elseif (0 < Y(k) && Y(k) <= 2*A)
        est_X(k) = 1*A;
    elseif (Y(k) > 2*A)
        est_X(k) = 3*A;
    end

end
