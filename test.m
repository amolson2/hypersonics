meshfile = 'fixed_left_tail_mesh_unshared.stl';

incoming_flow = [1 0 0];

rho_inf = [0.00104];
p_inf = [80.359];
alt = [49942];
v_inf = [2333.8];
M_inf = [7.076];
aoa = deg2rad([6.83]);
gamma = 1.4;

[pressures, cps] = supersonic_tangent_wedgie(meshfile, gamma, rho_inf,p_inf,v_inf,M_inf,aoa,incoming_flow);

figure(1)
c = pressures(:,3);
patch('Faces', mesh.Faces, 'Vertices', mesh.Vertices, 'CData', c, 'FaceColor', 'flat', 'EdgeColor', 'none')
view(3)
axis vis3d
axis equal
title('Fuselage Pressure Distribution')
cbar = colorbar;
cbar.Label.String = 'Pressure (Pa)';