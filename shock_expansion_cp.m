meshfile = 'fuselagemesh.stl';
mesh = readSurfaceMesh(meshfile);
M_inf = 7;
p_inf = 80;
rho_inf = 0.00104;
theta_nose = 0.5;
gamma = 1.4

%function [c_p, pressures] = shock_expansion_cp(mesh, gamma, M_inf, p_inf, rho_inf, theta_nose)

K_nose = M_inf * theta_nose;
beta_nose = theta_nose * ((gamma+1)/4 + sqrt(((gamma+1)/4)^2 + 1 / K_nose^2));
pnose_pinf = 1 + 2 * gamma / (gamma+1) * (((gamma + 1)/4 * K_nose + sqrt((gamma+1)/4)^2 * K_nose^2 + 1)^2 - 1);
p_nose = pnose_pinf * p_inf;
a_inf = sqrt(gamma * p_inf / rho_inf);
v_inf = M_inf * a_inf;
Mn_nose = sqrt(((gamma-1)*M_inf^2*sin(beta_nose)^2+2)/(2*gamma*M_inf^2*sin(beta_nose)^2-(gamma-1)));
M_nose = Mn_nose / sin(beta_nose - theta_nose);

[areas, planform, centers, thetas] = calculate_areas(mesh);

p_ratio = zeros(mesh.NumFaces, 1);
pressures = zeros(mesh.NumFaces, 1);
c_p = zeros(mesh.NumFaces, 1);

length(thetas)

for i=1:length(thetas)
    p_ratio(i) = expansion(M_nose, thetas(i) - theta_nose);
    pressures(i) = p_ratio(i) * p_nose;
    c_p(i) = (pressures(i) - p_inf) / (1 / 2 * rho_inf * v_inf^2);

    if i > 18500
        i
    end
end

%end