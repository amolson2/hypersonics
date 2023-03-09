function [c_p, pressures] = shock_expansion_cp(mesh, M_inf, p_inf, rho_inf, aoa, incoming_flow, theta_nose)

K_nose = M_inf * theta_nose;
beta_nose = theta_nose * ((gamma+1)/4 + sqrt(((gamma+1)/4)^2 + 1 / K_nose^2);
pnose_pinf = 1 + 2 * gamma / (gamma+1) * (((gamma + 1)/4 * K_nose + sqrt((gamma+1)/4)^2 * K_nose^2 + 1)^2 - 1);
p_nose = pnose_pinf * p_inf;
a_inf = sqrt(gamma * p_inf / rho_inf);
v_inf = M_inf * a_inf;

[areas, planform, centers, thetas] = calculate_areas(mesh);

p_ratio = zeros(mesh.NumFaces, 1);
rho_ratio = zeros(mesh.NumFaces, 1);
T_ratio = zeros(mesh.NumFaces, 1);
Machs = zeros(mesh.NumFaces, 1);
pressures = zeros(mesh.NumFaces, 1);
c_p = zeros(mesh.NumFaces, 1);

for i=1:length(thetas)
    [p_ratio(i), rho_ratio(i), T_ratio(i), Machs(i)] = expansion(M_nose, theta(i) - theta_nose);
    pressures(i) = p_ratio * p_nose;
    c_p(i) = (pressures(i) - p_inf) / (1 / 2 * rho_inf * v_inf^2);
end

end