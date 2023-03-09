function [pFactor, rhoFactor, Tfactor, M, beta] = oblique(M_inf, t)
    %Calculate conditions after a generic oblique shock
    gam = 1.4;

    beta = beta_calc(M_inf,t*pi/180,gam,0)*180/pi;

    Mn = sqrt(((gam-1)*M_inf^2*sind(beta)^2+2)/(2*gam*M_inf^2*sind(beta)^2-(gam-1)));
    M = Mn / sind(beta - t);
    
    pFactor =  (1 + (2*gam/(gam+1)) * (M_inf^2 * sind(beta)^2 - 1) );
    rhoFactor = (gam+1)*(M_inf^2*sind(beta)^2)/(2+(gam-1)*M_inf^2*sind(beta)^2);
    Tfactor = pFactor / rhoFactor;

end