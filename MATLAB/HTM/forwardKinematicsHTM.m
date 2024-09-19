function H = forwardKinematicsHTM(DH, frame)

%% Esta función calcula la cinemática directa hasta el n-ésimo marco. Devuelve una matriz de transformación homogénea de 4x4.
%{
DH: parámetros de Denavit Hartenberg como un ARREGLO (simbólico o numérico)
frame: el marco en movimiento, como un NÚMERO ENTERO, que se va a estudiar (considerado como la n-ésima fila del arreglo mencionado anteriormente)
%}

    % Inicialización de la matriz para almacenar el resultado del producto entre matrices
    H = eye(4);

    % Itera para alcanzar el marco deseado que se va a estudiar
    for i = 1 : frame

        % Ecuación de cinemática directa
        H = H * (Rz(DH(i, 1)) * Tz(DH(i, 2)) * Tx(DH(i, 3)) * Rx(DH(i, 4)));

    end
end