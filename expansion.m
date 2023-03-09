function pFactor = expansion(M_inf, t)
    %calculate conditions after a generic Prandtl-Meyer expansion
    gam = 1.4;
    Nu1 = sqrt((gam+1)/(gam-1))*atan(sqrt((gam-1)/(gam+1)*(M_inf^2 - 1))) - atan(sqrt(M_inf^2 - 1));
    M = invPrandtlMeyer(Nu1, t, gam);

    temp = (1+ (gam-1)/2*M_inf^2)/(1+ (gam-1)/2*M^2);
    pFactor = temp^(gam/(gam-1));
    %rhoFactor = temp^(1/(gam-1));
    %Tfactor = temp;
end
