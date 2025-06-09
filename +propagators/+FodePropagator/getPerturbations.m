function perturbations = getPerturbations(t, state, const, body)
    % Initialize perturbation in RTN frame
    perturbations = zeros(3, 1);

    % Add J2 perturbations
    J2Perturbations = propagators.FodePropagator.getJ2Perturbations(t, state, const, body);
  
    perturbations = perturbations + J2Perturbations;
end
