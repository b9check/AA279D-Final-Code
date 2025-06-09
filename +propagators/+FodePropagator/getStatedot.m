function statedot = getStatedot(t, state, const, body, perturbated)
    % Constant
    mu = const.(body).mu;
    
    % Extract state variables
    x = state(1);
    y = state(2);
    z = state(3);
    vx = state(4);
    vy = state(5);
    vz = state(6);

    % Precompute common terms
    r_vec = [x; y; z];
    v_vec = [vx; vy; vz];
    r = norm(r_vec);

    % Allocate output
    statedot = zeros(6, 1);

    % Two-body (Keplerian) acceleration
    accel_kepler = -mu / r^3 * r_vec;

    % Perturbations
    if nargin < 5 || perturbated
        perturbations = propagators.FodePropagator.getPerturbations(t, state, const, body);
    else
        perturbations = zeros(3, 1);
    end    

    % State derivative
    statedot(1:3) = v_vec;                    
    statedot(4:6) = accel_kepler + perturbations;
end
