function J2perturbations = getJ2Perturbations(t, state, const, body)
    % Constants
    mu = const.(body).mu; % Body gravitational parameter [m^3/s^2]
    R = const.(body).R;   % Body radius [m]
    J2 = const.(body).J2; % Body J2 constant [-]
    
    % Extract state variables
    a = state(1); % semi-major axis [m]
    e = state(2); % eccentricity [-]
    i = state(3); % inclination [rad]
    W = state(4); % RAAN (Omega) [rad]
    w = state(5); % argument of perigee [rad]
    f = state(6); % true anomaly [rad]

    % Precompute common terms
    r = a * (1 - e^2) / (1 + e * cos(f));
    sin_i = sin(i);
    cos_i = cos(i);
    sin_u = sin(w + f);
    cos_u = cos(w + f);
    
    % J2 factor
    J2_factor = -1.5 * mu * J2 * R^2 / r^4;

    % J2 perturbations in RTN frame
    f_R = J2_factor * (1 - 3 * sin_i^2 * sin_u^2);
    f_T = J2_factor * sin_i^2 * sin_u * cos_u;
    f_N = J2_factor * sin_i * cos_i * sin_u;
    
    % Output
    J2perturbations = [f_R; f_T; f_N];
end