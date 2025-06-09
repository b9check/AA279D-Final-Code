function statedot = getStatedot(t, state, const, body, perturbated)
    % Constant
    mu = const.(body).mu;

    % Extract state variables
    a = state(1); % semi-major axis [m]
    e = state(2); % eccentricity [-]
    i = state(3); % inclination [rad]
    W = state(4); % RAAN (Omega) [rad]
    w = state(5); % argument of perigee [rad]
    f = state(6); % true anomaly [rad]

    % Precompute common terms    
    p = a * (1 - e^2);
    r = p / (1 + e * cos(f));
    u = w + f;
    n = sqrt(mu / a^3);
    eta = sqrt(1-e^2);
    
    % Get perturbations
    if nargin < 5 || perturbated
        perturbations = propagators.KeplerianPropagator.getPerturbations(t, state, const, body);
    else
        perturbations = zeros(3, 1);
    end
    f_R = perturbations(1);
    f_T = perturbations(2);
    f_N = perturbations(3);
    
    % Allocate output
    statedot = zeros(6, 1);

    % Gauss' variational equations

    % da/dt
    statedot(1) = ((2*e*sin(f))/(n*eta))*f_R + ((2*a*eta)/(n*r))*f_T;

    % de/dt
    statedot(2) = ((eta*sin(f))/(n*a))*f_R + ((eta)/(n*a^2*e))*(((a^2*eta^2)/(r))-r)*f_T;

    % di/dt
    statedot(3) = ((r*cos(u))/(n*a^2*eta))*f_N;

    % dOmega/dt
    statedot(4) = ((r*sin(u))/(n*a^2*eta*sin(i)))*f_N;

    % domega/dt
    statedot(5) = ((-eta*cos(f))/(n*a*e))*f_R + ((eta)/(n*a*e))*((2+e*cos(f))/(1+e*cos(f)))*(sin(f))*f_T + ((-r*(1/tan(i))*sin(u))/(n*a^2*eta))*f_N;

    % df/dt
    statedot(6) = (sqrt(mu * p) / r^2) + (eta / (e * a * n)) * (cos(f) * f_R - (1 + r / p) * sin(f) * f_T);
end
