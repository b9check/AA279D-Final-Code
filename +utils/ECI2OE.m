function state_OE = ECI2OE(state_ECI, mu)    
    r_ECI = state_ECI(1:3);
    v_ECI = state_ECI(4:6);

    h_vec = cross(r_ECI, v_ECI);
    h = norm(h_vec);

    e_vec = cross(v_ECI, h_vec) / mu - r_ECI / norm(r_ECI);
    e = norm(e_vec);

    % True anomaly
    argf = dot(e_vec, r_ECI) / (e * norm(r_ECI));
    argf = min(max(argf, -1), 1);
    f = acos(argf);
    if dot(r_ECI, v_ECI) < 0
        f = 2*pi - f;
    end

    % Inclination
    cosi = h_vec(3) / h;
    cosi = min(max(cosi, -1), 1);
    i = acos(cosi);

    % Node vector
    N_vec = cross([0 0 1], h_vec);
    N = norm(N_vec);

    % RAAN
    if N ~= 0
        argW = N_vec(1) / N;
        argW = min(max(argW, -1), 1);
        W = acos(argW);
        if N_vec(2) < 0
            W = 2*pi - W;
        end
    else
        W = 0;
    end

    % Argument of periapsis
    if N ~= 0 && e > 1e-10
        argw = dot(N_vec, e_vec) / (N * e);
        argw = min(max(argw, -1), 1);
        w = acos(argw);
        if e_vec(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end

    % Semi-major axis
    a = 1 / (2 / norm(r_ECI) - norm(v_ECI)^2 / mu);
    state_OE = [a, e, i, W, w, f];
end
