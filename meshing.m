close all
clc

mesh = readSurfaceMesh('fuselagemesh.stl');
surfaceMeshShow(mesh)

%%

xC = zeros(mesh.NumFaces,1);
yC = zeros(mesh.NumFaces,1);
zC = zeros(mesh.NumFaces,1);

for i=1:mesh.NumFaces
    coords1 = mesh.Vertices(mesh.Faces(i,1), :);
    coords2 = mesh.Vertices(mesh.Faces(i,2), :);
    coords3 = mesh.Vertices(mesh.Faces(i,3), :);

    xC(i) = (coords1(1) + coords2(1) + coords3(1)) / 3;
    yC(i) = (coords1(2) + coords2(2) + coords3(2)) / 3;
    zC(i) = (coords1(3) + coords2(3) + coords3(3)) / 3;
end

%%
figure(1)
quiver3(xC, yC, zC, mesh.FaceNormals(:,1), mesh.FaceNormals(:,2), mesh.FaceNormals(:,3))
hold on
quiver3(0, 0, 0, 5, 0, 0)
hold off
