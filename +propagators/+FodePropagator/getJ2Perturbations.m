function J2perturbations = getJ2Perturbations(t, state, mu, R, J2)    
    % Extract state variables
    x = state(1);
    y = state(2);
    z = state(3);


    % Precompute common terms
    r_vec = [x; y; z];
    r = norm(r_vec);
    r2 = r^2;
    z2 = z^2;
    
    % J2 factor
    J2_factor = 1.5 * J2 * mu * R^2 / r^5;


    % J2 perturbations in ECI frame
    ax_J2 = J2_factor * x * (5 * z2 / r2 - 1);
    ay_J2 = J2_factor * y * (5 * z2 / r2 - 1);
    az_J2 = J2_factor * z * (5 * z2 / r2 - 3);
    
    % Output
    J2perturbations = [ax_J2; ay_J2; az_J2];
end