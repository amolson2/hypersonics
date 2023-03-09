
function [c_ptw, p] = supersonic_tangent_wedgie(mesh, giant_matrix)

p_inf = giant_matrix(2,:);
M_inf = giant_matrix(5,:);
aoa = giant_matrix(6,:);

incoming_flow = [1, 0, 0];
gamma = 1.4;
len = length(p_inf);

%% 
phi = zeros(mesh.NumFaces, len);
c_ptw = zeros(mesh.NumFaces, len);
p = zeros(mesh.NumFaces, len);
theta = zeros(mesh.NumFaces, len);

for j=1:len
    for i=1:mesh.NumFaces
        phi(i,j) = aoa(j) + pi - acos(dot(incoming_flow, mesh.FaceNormals(i, :)));
        theta(i,j) = pi/2 - phi(i,j);
       
        if theta(i,j) >= 0

            % Tangent Wedge
            [pFactor, rhoFactor, Tfactor, M,beta] = oblique(M_inf(j), theta(i,j)*180/pi);
            
            if isreal(pFactor) == false
                c_ptw(i,j) = 0;
                p(i,j) = p_inf(j);
            else
                p(i,j) = p_inf(j)*pFactor;
                c_ptw(i,j) = 2/gamma/M_inf(j)^2*(pFactor-1);
            end
        else 
            c_ptw(i,j) = 0;
            p(i,j) = p_inf(j);
        end

    end

end

