function [coefficients] = aerodynamic_coefficients(mesh, p_array, flow_conditions)
    n = mesh.FaceNormals;
    [S, planform, centers, thetas] = calculate_areas(mesh);
    % x = -centers;
    coefficients = zeros(5, 12);
    for j=1:12
        p = p_array(:,j);
        aoa = flow_conditions(6,j);
        rho_inf = flow_conditions(1,j);
        v_inf = flow_conditions(4,j);
        q = 0.5*rho_inf*v_inf^2;
        N = 0;
        A = 0;
        M = 0;
        for i=1:mesh.NumFaces
            N = N + dot(-p(i)*n(i,:)*S(i), [1 0 0]);
            A = A + dot(-p(i)*n(i,:)*S(i), [0 1 0]);
            L = N*cos(aoa) - A*sin(aoa);
            D = N*sin(aoa) + A*cos(aoa);
            % M = M + cross(-p(i)*n(i,:), x(i,:))*S(i);
        end
        coefficients(1,j) = N/(q*planform);
        coefficients(2,j) = A/(q*planform);
        % coefficients(3,j) = norm(M)/(q*planform^2);
        coefficients(4,j) = L/(q*planform);
        coefficients(5,j) = D/(q*planform);
    end
end