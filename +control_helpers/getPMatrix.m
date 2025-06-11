function P_bar = getPMatrix(chief_OE, current_QNSROE, desired_QNSROE, k, N)
    % Extract chief orbital elements
    a = chief_OE(1); % Semi-major axis (m)
    e = chief_OE(2); % Eccentricity (-)
    i = chief_OE(3); % Inclination (rad)
    W = chief_OE(4); % RAAN (rad)
    w = chief_OE(5); % Argument of perigee (rad)
    f = chief_OE(6); % True anomaly (rad)

    % Compute mean anomaly (M) from true anomaly (f)
    E = 2 * atan2(sqrt(1 - e) * sin(f/2), sqrt(1 + e) * cos(f/2)); % Eccentric anomaly
    M = E - e * sin(E); % Mean anomaly

    % Mean argument of latitude
    phi = w + M;

    % Compute tracking errors
    delta_rQNSOE = current_QNSROE - desired_QNSROE;
    delta_rQNSOE_ex = delta_rQNSOE(3);
    delta_rQNSOE_ey = delta_rQNSOE(4);
    delta_rQNSOE_ix = delta_rQNSOE(5);
    delta_rQNSOE_iy = delta_rQNSOE(6);

    % Compute fuel-optimal angles
    phi_ip = atan2(delta_rQNSOE_ey, delta_rQNSOE_ex); % In-plane optimal angle
    phi_oop = atan2(delta_rQNSOE_iy, delta_rQNSOE_ix); % Out-of-plane optimal angle
    phi_lambda = w; % Radial thrust at perigee (approximate)

    % Angular offsets
    J = phi - phi_ip;
    H = phi - phi_oop;
    K = phi - phi_lambda;

    % Ensure N is even and positive
    if mod(N, 2) ~= 0 || N < 2
        error('N must be an even positive integer >= 2');
    end

    % Construct P_bar
    P_bar = (1/k) * diag([cos(J)^N, cos(K)^N, cos(J)^N, cos(J)^N, cos(H)^N, cos(H)^N]);
end