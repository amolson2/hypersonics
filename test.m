meshfile = 'fixed_left_tail_mesh_unshared.stl';
mesh = readSurfaceMesh(meshfile);

incoming_flow = [1 0 0];

rho_inf = [0.00104];
p_inf = [80.359];
alt = [49942];
v_inf = [2333.8];
M_inf = [7.076];
aoa = deg2rad([6.83]);
gamma = 1.4;

[pressures, cps] = supersonic_tangent_wedgie(mesh, gamma, rho_inf,p_inf,v_inf,M_inf,aoa,incoming_flow);

for i = 1:mesh.NumFaces
    for j = 1:length(rho_inf)
        if isreal(pressures(i,j)) == false
            pressures(i,j) = -10000;
        end
    end
end
figure(1)
c = pressures(:,1);
patch('Faces', mesh.Faces, 'Vertices', mesh.Vertices, 'CData', c, 'FaceColor', 'flat', 'EdgeColor', 'none')
view(3)
axis vis3d
axis equal
title('Fuselage Pressure Distribution')
cbar = colorbar;
cbar.Label.String = 'Pressure (Pa)';