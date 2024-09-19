function x = axisAngle(H)

% Converts an Homogeneous Transformation Matrix into a Position and Axis - Angle vector

    % Angle of rotation of the vector: Φ(t)
    phi = acos((trace(H(1 : 3, 1 : 3)) - 1)/2);

    % Axis of rotation: n
    n = (1/(2 * sin(phi))) * [H(3,2) - H(2,3)
                              H(1,3) - H(3,1)
                              H(2,1) - H(1,2)];

    % Position and Axis - Angle vector: x = [r_x r_y_ r_z Φ * n]^T
    x = [H(1 : 3, 4)
         phi * n];
end