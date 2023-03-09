
function [p, c_ptw] = supersonic_tangent_wedgie(mesh, giant_matrix)

rho_inf = giant_matrix(1,:);
p_inf = giant_matrix(2,:);
v_inf = giant_matrix(4,:);
M_inf = giant_matrix(5,:);
aoa = giant_matrix(6,:);

incoming_flow = [1, 0, 0];
gamma = 1.4;
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
        phi(i,j) = aoa(j) + pi - acos(dot(incoming_flow, mesh.FaceNormals(i, :)));
        theta(i,j) = pi/2 - phi(i,j);
       
        if theta(i,j) >= 0

            % Tangent Wedge (from lecture slides)
            %K(i,j) = M_inf(j)*theta(i,j);
            %Ks(i,j) = (gamma + 1)/4*K(i,j) + sqrt(((gamma+1)/4)^2*K(i,j)^2+1);
            %c_ptw(i,j) = 4/(gamma + 1) * theta(i,j)^2*(Ks(i,j)^2-1)/K(i,j)^2;
            
            % Tangent Wedge (from our own intuition)
            [pFactor, rhoFactor, Tfactor, M,beta] = oblique(M_inf(j), theta(i,j)*180/pi);
            p(i,j) = p_inf(j)*pFactor;
           
            c_ptw(i,j) = (p(i,j) - p_inf(j))/(0.5*rho_inf(j)*v_inf(j)^2);
        else 
            c_ptw(i,j) = 0;
            p(i,j) = 0;
        end

%          if isreal(pFactor) == false
%              theta(i,j)*180/pi
%          end
    end

end

