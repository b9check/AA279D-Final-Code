function QNSROE = ROE2QNSROE(chief_OE, deputy_SROE)
% ROE2QNSROE  Convert singular relative orbital elements to quasi-nonsingular ROEs
%   chief_OE      = [a0; e0; i0; RAAN0; argp0; nu0]   (all angles in radians)
%   deputy_SROE   = [da; de; di; dOmega; domega; dnu] (all angles in radians)

% Extract chief orbital elements
e0      = chief_OE(2);
i0      = chief_OE(3);
RAAN0   = chief_OE(4);
argp0   = chief_OE(5);
nu0     = chief_OE(6);

% Extract singular relative elements
da      = deputy_SROE(1);
de      = deputy_SROE(2);
di      = deputy_SROE(3);
dRAAN   = deputy_SROE(4);
domega  = deputy_SROE(5);
dnu     = deputy_SROE(6);

% Convert delta-nu to delta-mean anomaly (dM) (all in radians)
E0 = 2 * atan( sqrt((1 - e0) / (1 + e0)) * tan(nu0 / 2) );
Ed = 2 * atan( sqrt((1 - e0) / (1 + e0)) * tan((nu0 + dnu) / 2) );
M0 = E0 - e0 * sin(E0);
Md = Ed - e0 * sin(Ed);
dM = Md - M0;

% Compute QNSROE components (all trigs in radians)
dlam = dM + domega + dRAAN * cos(i0);
dex  = de * cos(argp0) - e0 * domega * sin(argp0);
dey  = de * sin(argp0) + e0 * domega * cos(argp0);
dix  = di;
diy  = dRAAN * sin(i0);

% Assemble output
QNSROE = [da; dlam; dex; dey; dix; diy];
end
