function c_p = shock_expansion_cp(mesh, M_inf, p_inf, rho_inf, aoa, incoming_flow, theta_nose)

K_nose = M_inf * theta_nose;
beta_nose = theta_nose * ((gamma+1)/4 + sqrt(((gamma+1)/4)^2 + 1 / K_nose^2);
pnose_pinf = 1 + 2 * gamma / (gamma+1) * (((gamma + 1)/4 * K_nose + sqrt((gamma+1)/4)^2 * K_nose^2 + 1)^2 - 1);


end