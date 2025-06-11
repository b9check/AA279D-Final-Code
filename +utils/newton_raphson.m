function E = newton_raphson(M, e, epsilon)
    if nargin < 3
        epsilon = 1e-10;
    end

    E = M;
    max_iter = 1e5;

    for i = 1:max_iter
        f_E = E - e * sin(E) - M;
        f_prime_E = 1 - e * cos(E);
        increment = f_E / f_prime_E;
        E = E - increment;

        if abs(increment) <= epsilon
            break;
        end
    end
end