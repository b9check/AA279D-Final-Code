function rQNSOE = ECI2rQNSOE(chief_ECI, deputy_ECI, mu)
    N = size(chief_ECI, 1);
    rQNSOE = zeros(N, 6);
    for k = 1:N
        % Extract ECI states
        chief_state = chief_ECI(k, :)';
        deputy_state = deputy_ECI(k, :)';
        % Convert to Keplerian OE
        oe_chief = utils.ECI2OE(chief_state, mu);
        oe_deputy = utils.ECI2OE(deputy_state, mu);
        % Compute relative QNS
        rQNSOE(k, :) = utils.OE2rQNSOE(oe_chief, oe_deputy);
    end
end