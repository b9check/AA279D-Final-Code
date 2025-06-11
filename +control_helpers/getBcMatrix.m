function Bc = getBcMatrix(chief_OE, mu, J2, R)
     % Extract orbital elements
    a = chief_OE(1); e = chief_OE(2); i = chief_OE(3);
    W = chief_OE(4); w = chief_OE(5); f = chief_OE(6);

    % Derived quantities
    eta = sqrt(1 - e^2);
    nc = sqrt(mu / a^3);
    gamma = (3/4) * J2 * R^2 * sqrt(mu);
    kappa = gamma / (a^(7/2) * eta^4);

    % Trig values
    sin_wf = sin(w + f);
    cos_wf = cos(w + f);
    tan_i = tan(i);
    cos_f = cos(f);
    sin_f = sin(f);

    % Eccentricity vector components
    ex = e * cos(w);
    ey = e * sin(w);

    % Common factors
    denom = 1 + e * cos_f;
    ec2pf = (2 + e * cos_f); % shorthand

    % Compute B matrix
    Bc = zeros(6,3);
    Bc(1,1) = 2 * e * sin_f / eta;
    Bc(1,2) = 2 * (1 + e * cos_f) / eta;
    Bc(2,1) = -2 * eta^2 / denom;

    Bc(3,1) = eta * sin_wf;
    Bc(3,2) = eta * (ec2pf * cos_wf + ex) / denom;
    Bc(3,3) = eta * ey * sin_wf / (tan_i * denom);

    Bc(4,1) = -eta * cos_wf;
    Bc(4,2) = eta * (ec2pf * sin_wf + ey) / denom;
    Bc(4,3) = -eta * ex * sin_wf / (tan_i * denom);

    Bc(5,3) = eta * cos_wf / denom;

    Bc(6,3) = eta * sin_wf / denom;

    % Normalize
    Bc = Bc / (a * nc);
end
