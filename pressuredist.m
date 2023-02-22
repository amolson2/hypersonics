close all
clc
%%
meshfile = 'fuselagemesh.stl';
mesh = readSurfaceMesh(meshfile);
incoming_flow = [1 0 0];

% constants - example profile taken first entry in slide 10 of precept 1 slides
% in theory we could do this for all times eventually
rho_inf = [0.00104 0.0056 0.00732 0.01147 0.01859 0.02267 0.02493 0.02727 0.02734 0.02197 0.03580 0.0836];
p_inf = [80.359 392.181 502.252 760.605 1208.457 1465.014 1606.457 1753.112 1757.515 1421.038 2283.671 5206.06];
alt = [49942 37716.86 35947.24 33059.99 29936.43 28652.17 28040.09 27461.55 27444.96 28854.96 25720.26 20384.69];
v_inf = [2333.8 2337.68 2325.39 2288.09 2187.9 2136.93 2112 2066.64 1928.32 1446.28 1034.27 591.17];
M_1 = [7.076 7.465 7.501 7.507 7.253 7.104 7.03 6.888 6.427 4.806 3.461 2.002];

gamma = 1.4;
%% 
% simple way of calculating pressure distribution
phi = zeros(mesh.NumFaces, 1);
c_p = zeros(mesh.NumFaces, 1);
p = zeros(mesh.NumFaces, 12);

for j=1:12
    for i=1:mesh.NumFaces
        phi(i) = rad2deg(acos(dot(incoming_flow, mesh.FaceNormals(i, :))/norm(mesh.FaceNormals(i, :))));
        c_p(i) = 2*cos(phi(i))^2;
        p(i,j) = p_inf(j)+0.5*rho_inf(j)*v_inf(j)^2*c_p(i);
    end
end
%%
% higher fidelity - using shock relations from class
% this is still not right for some reason
theta = zeros(mesh.NumFaces, 1);
beta = zeros(mesh.NumFaces, 12);
c_p_m = zeros(mesh.NumFaces, 12);
p_m = zeros(mesh.NumFaces, 12);
for j=1:12
    for i=1:mesh.NumFaces
        theta(i) = 90 - phi(i); % is this right?
        beta(i,j) = beta_calc(M_1(j), theta(i), gamma, 0);
        c_p_m(i,j) = 4/(gamma+1)*(sin(beta(i,j))^2 - (1/M_1(j)^2));
        p_m(i,j) = p_inf(j)+0.5*rho_inf(j)*v_inf(j)^2*c_p_m(i,j);
    end
end