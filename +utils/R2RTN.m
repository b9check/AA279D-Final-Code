function R_ECI2RTN = R2RTN(r_sat_ECI, v_sat_ECI)
    % Compute the radial unit vector (R-hat)
    r_hat = r_sat_ECI / norm(r_sat_ECI);

    % Compute the angular momentum vector and then the normal unit vector (N-hat)
    h_vec = cross(r_sat_ECI, v_sat_ECI);
    n_hat = h_vec / norm(h_vec);

    % Compute the transverse unit vector (T-hat) via the cross product of N-hat and R-hat
    t_hat = cross(n_hat, r_hat);

    % Form the transformation (rotation) matrix from ECI to RTN
    R_ECI2RTN = [r_hat.'; t_hat.'; n_hat.'];
end