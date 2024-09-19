function M = cartesianMassMatrixHTM(D, J_g)
    
%% Esta función calcula la matriz de masa del robot en el espacio euclidiano
%{
    D: matriz de masa en el espacio de las articulaciones.
    J_g: matriz jacobiana geométrica.
%}

    % Inversa de la matriz jacobiana geométrica
    J_inv = inv(J_g);

    % Transformación de matrices
    M = transpose(J_inv) * D * J_inv;
    
end