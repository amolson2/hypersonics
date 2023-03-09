function [c_p, pressures] = shock_expansion_cp(mesh, giant_matrix)

rho_inff = giant_matrix(1,:);
p_inff = giant_matrix(2,:);
v_inff = giant_matrix(4,:);
M_inff = giant_matrix(5,:);
aoa = giant_matrix(6,:);

gamma = 1.4;
[areas, planform, centers, thetas] = calculate_areas(mesh);

pressures = zeros(length(rho_inff), mesh.NumFaces);
c_p = zeros(length(rho_inff), mesh.NumFaces);

for j = 1:length(rho_inff)
    rho_inf = rho_inff(j);
    p_inf = p_inff(j);
    v_inf = v_inff(j);
    M_inf = M_inff(j);

    theta_nose = max(thetas) - aoa(j);

    K_nose = M_inf * theta_nose;
    beta_nose = theta_nose * ((gamma+1)/4 + sqrt(((gamma+1)/4)^2 + 1 / K_nose^2));
    pnose_pinf = 1 + 2 * gamma / (gamma+1) * (((gamma + 1)/4 * K_nose + sqrt((gamma+1)/4)^2 * K_nose^2 + 1)^2 - 1);
    p_nose = pnose_pinf * p_inf;

    Mn_nose = sqrt(((gamma-1)*M_inf^2*sin(beta_nose)^2+2)/(2*gamma*M_inf^2*sin(beta_nose)^2-(gamma-1)));
    M_nose = Mn_nose / sin(beta_nose - theta_nose);
    
    p_ratio = zeros(mesh.NumFaces, 1);
    
    for i=1:length(thetas)
        p_ratio(i) = expansion(M_nose, theta_nose - thetas(i));
        pressures(j,i) = p_ratio(i) * p_nose;
        c_p(j,i) = (pressures(j,i) - p_inf) / (1 / 2 * rho_inf * v_inf^2);
    end

end

end
