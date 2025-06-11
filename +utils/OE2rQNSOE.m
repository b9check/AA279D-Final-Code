function rQNSOE = OE2rQNSOE(oe_chief, oe_deputy)
    % Extract chief elements
    a_0 = oe_chief(1);
    e_0 = oe_chief(2);
    i_0 = oe_chief(3);
    W_0 = oe_chief(4);
    w_0 = oe_chief(5);
    f_0 = oe_chief(6);
    M_0 = utils.TrueToMeanAnomaly(f_0, e_0);
    % Extract deputy elements
    a_1 = oe_deputy(1);
    e_1 = oe_deputy(2);
    i_1 = oe_deputy(3);
    W_1 = oe_deputy(4);
    w_1 = oe_deputy(5);
    f_1 = oe_deputy(6);
    M_1 = utils.TrueToMeanAnomaly(f_1, e_1);

    % Compute QNS-ROE
    delta_a  = (a_1 - a_0) / a_0;
    delta_lambda = (M_1 + w_1) - (M_0 + w_0) + (W_1 - W_0) * cos(i_0);
    delta_ex = e_1 * cos(w_1) - e_0 * cos(w_0);
    delta_ey = e_1 * sin(w_1) - e_0 * sin(w_0);
    delta_ix = i_1 - i_0;
    delta_iy = (W_1 - W_0) * sin(i_0);
    rQNSOE = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];
end