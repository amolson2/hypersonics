function [area_arr, planform, centers, thetas] = calculate_areas(mesh)

u_vec = [1, 0, 0];

xC = zeros(mesh.NumFaces,1);
yC = zeros(mesh.NumFaces,1);
zC = zeros(mesh.NumFaces,1);
areas = zeros(mesh.NumFaces,1);
theta = zeros(mesh.NumFaces,1);

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


    if mesh.FaceNormals(i, 3) > 0
        planform_area = planform_area + areas(i) * mesh.FaceNormals(i,3);
    end

end

area_arr = areas;
planform = planform_area;
centers = center_coord;
thetas = theta;

end