function plotOrbit3D_ECI(r, const, body)
    figure;
    hold on;
    axis equal;
    grid on;

    % Get body radius
    R_body = const.(body).R;

    % Plot the central body (sphere)
    [X, Y, Z] = sphere(100);
    surf(R_body * X / 1e3, R_body * Y / 1e3, R_body * Z / 1e3, ...
        'FaceColor', [0.8 0.8 0.8], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.7);

    % Plot orbit trajectory
    plot3(r(:,1)/1e3, r(:,2)/1e3, r(:,3)/1e3, 'b', 'LineWidth', 1.5);

    % Plot start and end points
    plot3(r(1,1)/1e3, r(1,2)/1e3, r(1,3)/1e3, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
    plot3(r(end,1)/1e3, r(end,2)/1e3, r(end,3)/1e3, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % End

    % Labels and title
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    title(['Orbit around ', upper(body), ' in ECI Frame']);

    % View settings
    view(3);
    camlight;
    lighting gouraud;

    legend(upper(body), 'Orbit Path', 'Start Point', 'End Point');
end
