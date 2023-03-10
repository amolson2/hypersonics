function [c_p_mod,c_p, p_mod] = corrected_newtonian(mesh, flight_cond, gamma)
    phi = zeros(mesh.NumFaces, 12);
    c_p = zeros(mesh.NumFaces, 12);
    p = zeros(mesh.NumFaces, 12);
    theta = zeros(mesh.NumFaces, 12);
    v = zeros(mesh.NumFaces, 12);

    incoming_flow = [1, 0, 0];

    for j=1:12
        rho_inf = flight_cond(1,j);
        p_inf = flight_cond(2,j);
        v_inf = flight_cond(4,j);
        M_inf = flight_cond(5,j);
        aoa = flight_cond(6,j);

        for i=1:mesh.NumFaces
            phi(i,j) = aoa + pi - acos(dot(incoming_flow, mesh.FaceNormals(i, :)));
            theta(i,j) = pi/2 - phi(i,j);
            c_p(i,j) = 2*cos(phi(i,j))^2;
            p(i,j) = p_inf+0.5*rho_inf*v_inf^2*c_p(i,j);
            v(i,j) = v_inf*cos(theta(i,j));
        end
    end

    % find stagnation pressure, c_p at stagnation, theta at stagnation
    pstag = zeros([1,12]);
    cpstag = zeros([1,12]);
    thetastag = zeros([1,12]);
    for i=1:12
        [min_phi, idx] = min(abs(v(:,i)));
        pstag(i) = p(idx,i);
        cpstag(i) = c_p(idx, i);
        thetastag(i) = theta(idx,i);
    end
    
    % Modified Newtonian
    c_p_mod = zeros(mesh.NumFaces, 12);
    p_mod = zeros(mesh.NumFaces, 12);
    for j=1:12
        rho_inf = flight_cond(1,j);
        p_inf = flight_cond(2,j);
        v_inf = flight_cond(4,j);
        M_inf = flight_cond(5,j);
        aoa = flight_cond(6,j);
        cpstag = (2/(gamma*M_inf^2))*(((gamma+1)^2*M_inf^2)/(4*gamma*M_inf^2 - (2*(gamma - 1))))^(gamma/(gamma-1))*((1 + (2*gamma)*(M_inf^2 - 1)/(gamma+1)) - 1);
        % cpstag = 2*thetastag(j)^2*((gamma+1)/4 + sqrt((gamma+1)/4) + (1/(M_inf*thetastag(j))^2));
        for i=1:mesh.NumFaces
            c_p_mod(i,j) = cpstag*(sin(theta(i,j)^2)/(sin(thetastag(j)^2)));
            p_mod(i,j) = p_inf+0.5*rho_inf*v_inf^2*c_p_mod(i,j);
        end
    end
end
