function mean_dot = getMeanStatedot(t, state, mu, R, J2, J2_perturbed)
    % state = [a; e; i; Omega; omega; nu]
    a     = state(1);
    e     = state(2);
    i     = state(3);
    Omega = state(4);
    omega = state(5);
    nu    = state(6);

    % Convert true anomaly ν → eccentric anomaly E
    E = 2 * atan( sqrt((1 - e) / (1 + e)) * tan(nu / 2) );

    % Convert eccentric anomaly E → mean anomaly M
    M = E - e * sin(E);

    % Common parameters
    p     = a * (1 - e^2);
    n     = sqrt(mu / a^3);
    cos_i = cos(i);

    if J2_perturbed
        % J2 secular perturbations (mean rates)
        dOmega_dt = -1.5 * n * J2 * (R^2 / p^2) * cos_i;
        domega_dt =  0.75 * n * J2 * (R^2 / p^2) * (5 * cos_i^2 - 1);
        dM_dt     =  n + 0.75 * n * J2 * (R^2 / p^2) * sqrt(1 - e^2) * (3 * cos_i^2 - 1);
    else
        % Two-body motion only
        dOmega_dt = 0;
        domega_dt = 0;
        dM_dt     = n;
    end

    % Time derivative of ν (based on chain rule through M)
    dE_dM   = 1 / (1 - e * cos(E));
    dnu_dE  = sqrt(1 - e^2) / (1 - e * cos(E));
    dnu_dt  = dnu_dE * dE_dM * dM_dt;

    % Mean state time derivative
    mean_dot = [
        0;           % da/dt (constant semi-major axis under J2)
        0;           % de/dt (no secular change under J2)
        0;           % di/dt (no secular change under J2)
        dOmega_dt;
        domega_dt;
        dnu_dt
    ];
end
