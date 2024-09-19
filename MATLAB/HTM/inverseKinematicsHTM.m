function theta = inverseKinematicsHTM(theta_0, L, r_d, K, frames, dt)

% This function calculates the inverse kinematics to the n - th frame. It returns an nx1 vector with the values of the joints to reach a desired pose H_d
%{
    theta_0: initial conditions of the joints
    L: length of the rigid bodies
    r_d: desired position as a VECTOR
    K: gain matrix (positive definite)
    frames: the moving frame, as an INTEGER NUMBER, to be studied (considered as the n - th row of the Denavit Hartenberg array)
    dt: step for the integration algorithm
%}

    % Start number of iterations
    i = 1;

    % Set initial conditions
    theta(:, i) = theta_0;

    % Iterates until the number of iterations is equal to 1500 (to limit the algorithm)
    while i <= 15000

        % Set the current values for the Denavit Hartenberg parameters
        DH = denavitHartenberg(theta(:, i), L);

        % Calculate forward kinematics
        H = forwardKinematicsHTM(DH, frames);
        
        % Calculate error functions (position + axis - angle vector)
        e = r_d - H(1 : 3, 4);

        % If the maximum absolute value of the error is less equal than 1 mm or 1 milirad
        if max(abs(e)) <= 1e-3
            
            % Show message
            fprintf("\nSOLUTION FOUND! Iterations: %d\n\n", i);

            % Stop the program
            break

        % Else,
        else
            
            % Calculate the Jacobian Matrix
            J = jacobianMatrixHTM(DH, theta(:, i), frames);

            % Calculate the inverse kinematics
            theta(:, i + 1) = theta(:, i) + (pinv(J) * K * e * dt);

        end

        % Counter plus one
        i = i + 1;
    end
end