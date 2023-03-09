function M2 = invPrandtlMeyer(nu1, theta, gamma)
    % using bisection algorithm to numerically solve for nu inverse and
    % find M2
    min = 0;
    max = 300;
    e = 2;
    nu = nu1 + theta;
    while abs(e) > 0.0001
        M2 = (min + max)/2;
        nu_guess = sqrt((gamma+1)/(gamma-1)) * atan(sqrt((gamma-1)/(gamma+1)*(M2^2-1))) - atan(sqrt(M2^2-1));
        e = nu_guess - nu;
        if e < 0
            min = M2;
        else
            max = M2;
        end
    end
end