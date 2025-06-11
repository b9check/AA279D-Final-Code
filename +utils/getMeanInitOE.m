function mean_OE0 = getMeanInitOE(osc_OE0, e, mu, Re, J2)
    % unpack
    a0     = osc_OE0(1);
    e0     = osc_OE0(2);
    i0     = osc_OE0(3);
    Omega0 = osc_OE0(4);
    omega0 = osc_OE0(5);
    nu0    = osc_OE0(6);

    % auxiliary
    p0   = a0*(1 - e0^2);
    eta0 = sqrt(1 - e0^2);
    q1   = e0*cos(omega0);
    q2   = e0*sin(omega0);
    theta0 = omega0 + nu0;
    E0   = 2*atan2( sqrt(1-e0)*tan(nu0/2), sqrt(1+e0) );
    M0   = E0 - e0*sin(E0);

    % common J2 factor
    K = J2*(Re^2)/(2*p0^2);

    % --- LP corrections (E.1–E.7) ---
    delta_a_lp = 0;
    % (E.2)
    delta_lambda_lp = q1*q2*sin(i0)^2/(8*a0^2*eta0^2*(1+eta0))*(1 - 5*cos(i0)^2) ...
                    + q1*q2/(16*a0^2*eta0^4)*(3 - 2.5*cos(i0)^2 - 14*cos(i0)^4 - 20*cos(i0)^6);
    % (E.3) not needed for e
    % (E.4) not needed for e
    % (E.5–E.6)
    delta_q1_lp = -q1*sin(i0)^2/(16*a0^2*eta0^2)*(1 - 5*cos(i0)^2) ...
                  - q1*q2^2/(16*a0^2*eta0^4)*(3 - 5*cos(i0)^2 - 28*cos(i0)^4 - 40*cos(i0)^6);
    delta_q2_lp =  q2*sin(i0)^2/(16*a0^2*eta0^2)*(1 - 5*cos(i0)^2) ...
                  + q1^2*q2/(16*a0^2*eta0^4)*(3 - 5*cos(i0)^2 - 28*cos(i0)^4 - 40*cos(i0)^6);
    % δe_lp from q‐vector linearization:
    delta_e_lp = (q1*delta_q1_lp + q2*delta_q2_lp)/e0;
    % (E.7) not needed for e

    % --- SP1 corrections (E.8–E.14): only δe_sp1 needed here
    delta_e_sp1 = K * e0 * (3*cos(i0)^2 - 1) * sin(2*E0) / 2;

    % --- SP2 corrections (E.15–E.21): none for e
    delta_e_sp2 = 0;

    % pack mean elements
    a_mean     = a0     - (delta_a_lp);
    e_mean     = e0     - (delta_e_lp + delta_e_sp1 + delta_e_sp2);
    i_mean     = i0;     % similar fix can be done for i if needed
    Omega_mean = Omega0; % ...
    omega_mean = omega0;
    M_mean     = M0;
    f_mean = utils.MeanToTrueAnomaly(M_mean, e);

    mean_OE0 = [a_mean; e_mean; i_mean; Omega_mean; omega_mean; f_mean];
end
