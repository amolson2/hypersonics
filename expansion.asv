function [pFactor, rhoFactor, Tfactor, M] = expansion(M_inf, t)
    %calculate conditions after a generic Prandtl-Meyer expansion
    gam = 1.4;
    t = t*pi/180;
    syms M_temp;
    Nu1 = sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1)*(M_inf^2 - 1))) - atan(sqrt(M_inf^2 - 1));
    eqnNu = sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1)*(M_temp^2 - 1))) - atan(sqrt(M_temp^2 - 1)) == Nu1 + t;
    M = vpasolve(eqnNu, M_temp,4);
    M = double(M);

    temp = (1+ (gam-1)/2*M_inf^2)/(1+ (gam-1)/2*M^2);
    pFactor = temp^(gam/(gam-1));
    rhoFactor = temp^(1/(gam-1));
    Tfactor = temp;
end
