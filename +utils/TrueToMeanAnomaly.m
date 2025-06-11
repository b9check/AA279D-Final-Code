function M = TrueToMeanAnomaly(f, e)
E = 2 * atan( sqrt((1 - e)/(1 + e)) * tan(f/2) );
E = mod(E, 2*pi);
M = mod(E - e*sin(E), 2*pi);
end
