function Phi_J2_qns = getQNS_J2_STM(OE, tau, mu, J2, R)
    % Extract orbital elements
    a = OE(1); e = OE(2); i = OE(3);
    W = OE(4); w = OE(5);

    % Derived quantities
    eta = sqrt(1 - e^2);
    kappa = (3/4) * J2 * R^2 * sqrt(mu) / (a^(7/2) * eta^4);

    % Auxiliary terms
    E = 1 + eta;
    F = 4 + 3 * eta;
    G = 1 / eta^2;
    P = 3 * cos(i)^2 - 1;
    Q = 5 * cos(i)^2 - 1;
    R = cos(i);
    S = sin(2*i);
    T = sin(i)^2;

    % Perifocal eccentricity vector components
    exi = e * cos(w);
    eyi = e * sin(w);

    % Secular rates
    w_dot = kappa * Q;
    
    % Final perigee argument
    omega_f = w + w_dot * tau;
    exf = e * cos(omega_f);
    eyf = e * sin(omega_f);

    % Trig terms
    cos_wt = cos(w_dot * tau);
    sin_wt = sin(w_dot * tau);

    % Construct STM
    Phi_J2_qns = zeros(6,6);
    Phi_J2_qns(1,1) = 1;
    Phi_J2_qns(2,1) = -((3/2)*sqrt(mu/a^3) + (7/2)*kappa*E*P) * tau;
    Phi_J2_qns(2,2) = 1;
    Phi_J2_qns(2,3) = kappa * exi * F * G * P * tau;
    Phi_J2_qns(2,4) = kappa * eyi * F * G * P * tau;
    Phi_J2_qns(2,5) = -kappa * F * S * tau;

    Phi_J2_qns(3,1) = (7/2) * kappa * eyf * Q * tau;
    Phi_J2_qns(3,3) = cos_wt - 4 * kappa * exi * eyf * G * Q * tau;
    Phi_J2_qns(3,4) = -sin_wt - 4 * kappa * eyi * eyf * G * Q * tau;
    Phi_J2_qns(3,5) = 5 * kappa * eyf * S * tau;

    Phi_J2_qns(4,1) = -(7/2) * kappa * exf * Q * tau;
    Phi_J2_qns(4,3) = sin_wt + 4 * kappa * exi * exf * G * Q * tau;
    Phi_J2_qns(4,4) = cos_wt + 4 * kappa * eyi * exf * G * Q * tau;
    Phi_J2_qns(4,5) = -5 * kappa * exf * S * tau;

    Phi_J2_qns(5,5) = 1;
    Phi_J2_qns(6,1) = (7/2) * kappa * S * tau;
    Phi_J2_qns(6,3) = -4 * kappa * exi * G * S * tau;
    Phi_J2_qns(6,4) = -4 * kappa * eyi * G * S * tau;
    Phi_J2_qns(6,5) = 2 * kappa * T * tau;
    Phi_J2_qns(6,6) = 1;
end