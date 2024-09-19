function C = centrifugalCoriolis(D, theta, dtheta)

%% Calculates the Centrifugal and Coriolis effects' matrix

%{
    D: mass matrix (SYMBOLIC)
    theta: vector of generalized coordinates (in RADIANS)
    dtheta: vector of generalized coordinates (in RADIANS / SECOND)
%}
    
    % Get number of columns of inertia matrix
    [~, n] = size(D);
    
    % Initialize matrix
    C = sym(zeros(n));

    % Iterate through all the columns of inertia matrix
    for k = 1 : n

        % Calculate derivative with respect to generalized coordinates
        dD_dtheta = jacobian(D(:, k), theta);

        % Calculate matrix
        C = C + ((dD_dtheta - (0.5 * transpose(dD_dtheta))) * dtheta(k));
    end

end
