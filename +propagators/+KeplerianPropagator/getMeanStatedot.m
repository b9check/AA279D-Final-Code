function mean_dot = getMeanStatedot(t, state, mu, R, J2)
    % state = [a; e; i; Omega; omega; nu]   (osculating a,e,i,Ω,ω,ν)
    a     = state(1);
    e     = state(2);
    i     = state(3);
    Omega = state(4);
    omega = state(5);
    nu    = state(6);

    % 1) convert ν -> eccentric anomaly E
    E = 2*atan( sqrt((1-e)/(1+e)) * tan(nu/2) );

    % 2) convert E -> mean anomaly M
    M = E - e*sin(E);

    % common parameters
    p     = a*(1 - e^2);
    n     = sqrt(mu/a^3);
    cos_i = cos(i);

    % secular J2‐averaged rates (in mean elements)
    da_dt     = 0;
    de_dt     = 0;
    di_dt     = 0;
    dOmega_dt = -1.5 * n * J2 * (R^2 / p^2) * cos_i;
    domega_dt =  0.75 * n * J2 * (R^2 / p^2) * (5*cos_i^2 - 1);
    dM_dt     =  n + 0.75 * n * J2 * (R^2 / p^2) * sqrt(1 - e^2) * (3*cos_i^2 - 1);

    % 3) convert dM/dt -> dν/dt via dν/dE * dE/dM * dM/dt:
    dE_dM   = 1 / (1 - e*cos(E));
    dnu_dE  = sqrt(1 - e^2) / (1 - e*cos(E));
    dnu_dt  = dnu_dE * dE_dM * dM_dt;

    % output derivatives of [a, e, i, Omega, omega, nu]
    mean_dot = [
        da_dt;
        de_dt;
        di_dt;
        dOmega_dt;
        domega_dt;
        dnu_dt
    ];
end
