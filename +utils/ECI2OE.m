function state_OE = ECI2OE(state_ECI, const, body)
    % Gravitational parameter
    mu = const.(body).mu;
    
    % Extrtact position and velocity
    r_ECI = state_ECI(1:3);
    v_ECI = state_ECI(4:6);
    
    % Specific angular momentum
    h_vec = cross(r_ECI, v_ECI);
    h = norm(h_vec);

    % Eccentricity vector and magnitude
    e_vec = cross(v_ECI, h_vec) / mu - r_ECI / norm(r_ECI);
    e = norm(e_vec);

    % True anomaly (f)
    arg = dot(e_vec, r_ECI) / (e * norm(r_ECI));
    arg = min(max(arg, -1), 1); % clamp to [-1, 1]
    f = acos(arg);
    if dot(r_ECI, v_ECI) < 0
        f = 2 * pi - f;
    end

    % Inclination (i)
    i = acos(h_vec(3) / h);

    % Node vector
    N_vec = cross([0, 0, 1], h_vec);
    N = norm(N_vec);

    % Right ascension of ascending node (W)
    W = acos(N_vec(1) / N);
    if N_vec(2) < 0
        W = 2 * pi - W;
    end

    % Argument of periapsis (w)
    w = acos(dot(N_vec, e_vec) / (N * e));
    if e_vec(3) < 0
        w = 2 * pi - w;
    end

    % Semi-major axis (a)
    a = 1 / (2 / norm(r_ECI) - norm(v_ECI)^2 / mu);
    state_OE = [a, e, i, W, w, f];
end
