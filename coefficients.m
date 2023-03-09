function [f_x, f_y, f_z, areas, forces]=coefficients(mesh, flight_conds)

    len = length(flight_conds(6,:));

    [cps, pressures] = supersonic_tangent_wedgie(mesh, flight_conds);
    
    [areas, planform, centers, thetas] = calculate_areas(mesh);

    norms = mesh.FaceNormals;
    
    forces = pressures .* areas;

    for j = 1:len
        f_x(j) = -1*sum(forces(:,j) .* norms(:,1));
        f_y(j) = -1*sum(forces(:,j) .* norms(:,2));
        f_z(j) = -1*sum(forces(:,j) .* norms(:,3));
    end

    %c_n = f_x / 

end

