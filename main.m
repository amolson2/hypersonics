close all
clc
%%
nose_file = 'fuselagemesh.stl';
nose_LE_file;
wing_LE_file;
wing_WS_file;
wing_LS_file;
fuselage_file = 'fuselagemesh.stl';
tail_LE_file;
tail_file;
flap_file;

nose = readSurfaceMesh(nosefile);
incoming_flow = [1 0 0];

rho_inf = [0.00104 0.0056 0.00732 0.01147 0.01859 0.02267 0.02493 0.02727 0.02734 0.02197 0.03580 0.0836];
p_inf = [80.359 392.181 502.252 760.605 1208.457 1465.014 1606.457 1753.112 1757.515 1421.038 2283.671 5206.06];
alt = [49942 37716.86 35947.24 33059.99 29936.43 28652.17 28040.09 27461.55 27444.96 28854.96 25720.26 20384.69];
v_inf = [2333.8 2337.68 2325.39 2288.09 2187.9 2136.93 2112 2066.64 1928.32 1446.28 1034.27 591.17];
M_inf = [7.076 7.465 7.501 7.507 7.253 7.104 7.03 6.888 6.427 4.806 3.461 2.002];
aoa = deg2rad([6.83 12 12 12 12 3.62 1.63 -0.66 -1.63 0.51 0.97 1.31]);
flight_conds = [rho_inf; p_inf; alt; v_inf; M_inf; aoa];

gamma = 1.4;
%% Newtonian method

%% Corrected Newtonian

%% Shock Expansion

%% Tangent Wedge

%% Make Figures
figure(1)
c = c_p(:,12);
patch('Faces', nose.Faces, 'Vertices', nose.Vertices, 'CData', c, 'FaceColor', 'flat', 'EdgeColor', 'none')
view(3)
axis vis3d
axis equal
title('Fuselage Pressure Distribution')
cbar = colorbar;
cbar.Label.String = 'Pressure (Pa)';

figure(2)
c = v(:,12);
patch('Faces', nose.Faces, 'Vertices', nose.Vertices, 'CData', c, 'FaceColor', 'flat', 'EdgeColor', 'none')
view(3)
axis vis3d
axis equal
title('Fuselage Velocity Distribution')
cbar = colorbar;
cbar.Label.String = 'Velocity (m/s)';