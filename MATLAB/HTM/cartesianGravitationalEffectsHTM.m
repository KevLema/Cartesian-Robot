function G_c = cartesianGravitationalEffectsHTM(G, J_g)
    
%% Esta función calcula los efectos gravitacionales del robot en el espacio euclidiano
%{
    G: vector of gravitational effects
    J_g: geometric Jacobian matrix
%}

    % Inverse of geometric Jacobian matrix
    J_inv = inv(J_g);

    % Matrix transformation
    G_c = transpose(J_inv) * G;

end