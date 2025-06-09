function x_RTN = ECI2RTN(chief_ECI, deputy_ECI)
r0 = chief_ECI(1:3)';
v0 = chief_ECI(4:6)';
r1 = deputy_ECI(1:3)';
v1 = deputy_ECI(4:6)';

% ECI TO CHIEF RTN ROTATION MATRIX
R0_hat = r0 / norm(r0);
N0 = cross(r0, v0);
N0_hat = N0 / norm(N0);
T0 = cross(N0_hat, R0_hat);
T0_hat = T0 / norm(T0);
R_ECI2RTN = [R0_hat'; T0_hat'; N0_hat'];

% Relative position in RTN frame
r_rel = r1 - r0;
rRTN = R_ECI2RTN * r_rel;

% Relative velocity
omega_RTN = cross(r0, v0) / norm(r0)^2;
v_rel_inertial = v1 - v0;
v_rel_rotating = v_rel_inertial - cross(omega_RTN, r_rel);
vRTN = R_ECI2RTN * v_rel_rotating;

x_RTN = [rRTN; vRTN];

end
