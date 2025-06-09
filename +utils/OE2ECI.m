function state_ECI = OE2ECI(state_OE,const, body)
    
    % Extract state variables
    i = state_OE(3); % inclination [rad]
    W = state_OE(4); % RAAN (Omega) [rad]
    w = state_OE(5); % argument of perigee [rad]

    % Get perifocal state
    state_perifocal = utils.OE2Perifocal(state_OE,const, body);

    % Rotation matrices
    R3_W = [ cos(W), -sin(W), 0;
             sin(W),  cos(W), 0;
                  0,       0, 1];

    R1_i = [1,      0,       0;
            0, cos(i), -sin(i);
            0, sin(i),  cos(i)];

    R3_w = [ cos(w), -sin(w), 0;
             sin(w),  cos(w), 0;
                  0,       0, 1];

    % Total rotation matrix 
    R_tot = R3_W * R1_i * R3_w;

    % Rotate position and velocity vectors
    r_ECI = R_tot * state_perifocal(1:3);
    v_ECI = R_tot * state_perifocal(4:6);

    % Output state vector
    state_ECI = [r_ECI; v_ECI];
end