function f = MeanToTrueAnomaly(M, e)
    % Newton-Raphson to solve for E
    E = M;
    tol = 1e-12;
    max_iter = 1e5;
    for i = 1:max_iter
        f_E = E - e * sin(E) - M;
        f_prime_E = 1 - e * cos(E);
        dE = f_E / f_prime_E;
        E = E - dE;
        if abs(dE) < tol
            break;
        end
    end

    % Convert E to true anomaly
    f = 2 * atan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2));
    f = mod(f, 2*pi);
end
