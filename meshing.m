clear
mesh = readSurfaceMesh('file.stl');
%surfaceMeshShow(mesh)

%%

mach = 6;
u_vec = [1, 0, 0];
v_vec = [0, 0, 1];
gamma = 1.4;

xC = zeros(mesh.NumFaces,1);
yC = zeros(mesh.NumFaces,1);
zC = zeros(mesh.NumFaces,1);
areas = zeros(mesh.NumFaces,1);
theta = zeros(mesh.NumFaces, 1);
beta = zeros(mesh.NumFaces, 1);

%%
wetted_area = 0;
non_wetted_area = 0;
planform_area = 0;
non_planform_area = 0;

for i=1:mesh.NumFaces
    coords1 = mesh.Vertices(mesh.Faces(i,1), :);
    coords2 = mesh.Vertices(mesh.Faces(i,2), :);
    coords3 = mesh.Vertices(mesh.Faces(i,3), :);

    xC(i) = (coords1(1) + coords2(1) + coords3(1)) / 3;
    yC(i) = (coords1(2) + coords2(2) + coords3(2)) / 3;
    zC(i) = (coords1(3) + coords2(3) + coords3(3)) / 3;

    AB = [coords2(1) - coords1(1), coords2(2) - coords1(2), coords2(3) - coords1(3)];
    AC = [coords3(1) - coords1(1), coords3(2) - coords1(2), coords3(3) - coords1(3)];

    center_coord = [xC, yC, zC];
    areas(i) = 0.5 * norm(cross(AB, AC));

    theta(i) = u_vec(1) * mesh.FaceNormals(i,1) + u_vec(2) * mesh.FaceNormals(i,2) + u_vec(3) * mesh.FaceNormals(i,3);
    if theta(i) < 0
        wetted_area = wetted_area + areas(i);
    else
        non_wetted_area = non_wetted_area + areas(i);
    end

    if mesh.FaceNormals(i, 3) > 0
        planform_area = planform_area + areas(i) * mesh.FaceNormals(i,3);
    else
        non_planform_area = non_planform_area + areas(i) * mesh.FaceNormals(i,3);
    end

    %syms bet
    %eqn = tan(bet-theta) / tan(bet) == (gamma - 1) / (gamma + 1);
    %beta(i) = solve(eqn,bet);
end

%%
quiver3(xC, yC, zC, mesh.FaceNormals(:,1), mesh.FaceNormals(:,2), mesh.FaceNormals(:,3))
axis equal

%%

beta = beta_calc(mach, theta, gamma, 0);
c_p = 4 / (gamma+1) * sin(beta).^2; 

%%
function Beta=beta_calc(M,theta,gamma,n)
    mu=asin(1/M);                   % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta);
    d=sqrt(4*(1-3*a.*b).^3./((27*a.^2*c+9*a.*b-2).^2)-1);
    Beta=atan((b+9*a*c)./(2*(1-3*a.*b))-(d.*(27*a.^2*c+9*a.*b-2))./(6*a.*(1-3*a.*b)).*tan(n*pi/3+1/3*atan(1./d)))*180/pi;
end