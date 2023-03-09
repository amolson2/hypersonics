
function [p, c_ptw] = supersonic_tangent_wedgie(meshfile, gamma, rho_inf, p_inf, v_inf, M_inf, aoa, incoming_flow)

mesh = readSurfaceMesh(meshfile);
len = length(rho_inf);

%% 
% Newtonian method
phi = zeros(mesh.NumFaces, len);
c_ptw = zeros(mesh.NumFaces, len);
p = zeros(mesh.NumFaces, len);
theta = zeros(mesh.NumFaces, len);
K = zeros(mesh.NumFaces, len);
Ks = zeros(mesh.NumFaces, len);

for j=1:len
    for i=1:mesh.NumFaces
        phi(i,j) = aoa(j) + acos(dot(incoming_flow, mesh.FaceNormals(i, :))/norm(mesh.FaceNormals(i, :)))+pi;
        theta(i,j) = pi/2 - phi(i,j);
        
        % Tangent Wedge (from lecture slides)
        K(i,j) = M_inf(j)*theta(i,j);
        Ks(i,j) = (gamma + 1)/4*K(i,j) + sqrt(((gamma+1)/4)^2*K(i,j)^2+1);
        c_ptw(i,j) = 4/(gamma + 1) * phi(i,j)^2*(Ks(i,j)^2-1)/K(i,j)^2;
        
        % Tangent Wedge (from our own intuition)
        %[pFactor, rhoFactor, Tfactor, M] = oblique(M_inf(j), theta(i,j)*180/pi);
        pFactor = 1;
        p(i,j) = p_inf(j)*pFactor;
        c_ptw(i,j) = (p(i,j) - p_inf(j))/(0.5*rho_inf(j)*v_inf(j)^2);

        if mod(i,100) == 0
            i
        end
    end

end
