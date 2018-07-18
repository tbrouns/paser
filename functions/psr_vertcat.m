function M = psr_vertcat(M1,M2)

    n = size(M1,2);
    m = size(M2,2);
    dn = n - m;

    if     (dn < 0); M1 = [M1,NaN(size(M1,1),abs(dn))];
    elseif (dn > 0), M2 = [M2,NaN(size(M2,1),abs(dn))];
    end

    M = [M1;M2];

end