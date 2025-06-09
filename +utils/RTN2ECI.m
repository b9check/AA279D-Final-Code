function [rECI, vECI] = RTN2ECI(r_chief_ECI, v_chief_ECI, r_deputy_RTN, v_deputy_RTN)
r0 = r_chief_ECI;
v0 = v_chief_ECI;
r1 = r_deputy_RTN;
v1 = v_deputy_RTN;

% CHIEF RTN TO ECI ROTATION MATRIX
R0_hat = r0 / norm(r0);
N0 = cross(r0, v0);
N0_hat = N0 / norm(N0);
T0 = cross(N0_hat, R0_hat);
T0_hat = T0 / norm(T0);
R_RTN2ECI = [R0_hat'; T0_hat'; N0_hat']';

% CALCULATE DEPUTY RTN COORDS
rECI = r0 + R_RTN2ECI * r1;
vECI = v0 + R_RTN2ECI * v1;

end
