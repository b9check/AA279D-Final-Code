function state_perifocal = OE2Perifocal(state_OE,const, body)
    % Constant
    mu = const.(body).mu;

    % Extract state variables
    a = state_OE(1); % semi-major axis [m]
    e = state_OE(2); % eccentricity [-]
    f = state_OE(6); % true anomaly [rad]
    
    
    % Compute the distance
    r = a * (1 - e^2) / (1 + e * cos(f));

    % Position in perifocal frame
    r_perifocal = [r * cos(f);
                   r * sin(f);
                   0];

    % Velocity in perifocal frame
    h = sqrt(mu * a * (1 - e^2)); % angular momentum
    v_perifocal = [-mu / h * sin(f);
                    mu / h * (e + cos(f));
                    0];

    % Output state vector
    state_perifocal = [r_perifocal; v_perifocal];
end
