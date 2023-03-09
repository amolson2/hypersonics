function [pFactor, rhoFactor, Tfactor, M] = oblique(M_inf, t)
    %Calculate conditions after a generic oblique shock
    gam = 1.4;
    syms b;
    
    eqnBeta = tand(b-t)/tand(b) == (2+(gam-1)*M_inf^2*sind(b)^2)/((gam+1)*M_inf^2*sind(b)^2);
    beta = vpasolve(eqnBeta, b, 20);
    beta = double(beta);

    Mn = sqrt(((gam-1)*M_inf^2*sind(beta)^2+2)/(2*gam*M_inf^2*sind(beta)^2-(gam-1)));
    M = Mn / sind(beta - t);
    
    pFactor =  (1 + (2*gam/(gam+1)) * (M_inf^2 * sind(beta)^2 - 1) );
    rhoFactor = (gam+1)*(M_inf^2*sind(beta)^2)/(2+(gam-1)*M_inf^2*sind(beta)^2);
    Tfactor = pFactor / rhoFactor;

end