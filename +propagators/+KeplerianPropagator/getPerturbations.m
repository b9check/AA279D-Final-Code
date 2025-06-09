function perturbations = getPerturbations(t, state, mu, R, J2)
    % Initialize perturbation in RTN frame
    perturbations = zeros(3, 1);

    % Add J2 perturbations
    J2Perturbations = propagators.KeplerianPropagator.getJ2Perturbations(t, state, mu, R, J2);
  
    perturbations = perturbations + J2Perturbations;
end
