close all
clc
%%
meshfile = 'fuselagemesh.stl';
mesh = readSurfaceMesh(meshfile);
incoming_flow = [1 0 0];

% constants - example profile taken first entry in slide 10 of precept 1 slides
% in theory we could do this for all times eventually
rho_inf = 0.00104;
p_inf = 80.359;
alt = 49942;
v_inf = 2333.8;
gamma = 1.4;
M_1 = 7.076;
%% 
% simple way of calculating pressure distribution
phi = zeros(mesh.NumFaces, 1);
c_p = zeros(mesh.NumFaces, 1);
p = zeros(mesh.NumFaces, 1);
for i=1:mesh.NumFaces
    phi(i) = rad2deg(acos(dot(incoming_flow, mesh.FaceNormals(i, :))/norm(mesh.FaceNormals(i, :))));
    c_p(i) = 2*cos(phi(i))^2;
    p(i) = p_inf+0.5*rho_inf*v_inf^2*c_p(i);
end
%%
% higher fidelity - using shock relations from class
% this is still not right for some reason
theta = zeros(mesh.NumFaces, 1);
beta = zeros(mesh.NumFaces, 1);
c_p_m = zeros(mesh.NumFaces, 1);
p_m = zeros(mesh.NumFaces, 1);
for i=1:mesh.NumFaces
    theta(i) = phi(i); % is this right?
    beta(i) = beta_calc(M_1, theta(i), gamma, 0);
    c_p_m(i) = 4/(gamma+1)*(sin(beta(i))^2 - (1/M_1^2));
    p_m(i) = p_inf+0.5*rho_inf*v_inf^2*c_p_m(i);
end

