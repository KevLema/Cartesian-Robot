function J = jacobianMatrixHTM(DH, theta, frame)

%% Calcula la Matriz Jacobiana Geométrica dada los parámetros de Denavit - Hartenberg y las posiciones de las articulaciones

%{

DH: Parámetros de Denavit-Hartenberg como un ARRAY (simbólico o numérico)
theta: vector de coordenadas generalizadas (en RADIANTES)
frames: el número de filas en tu array de Denavit - Hartenberg

%}

    % Obtener el número de articulaciones
    [n, ~] = size(theta);

    % Si el vector de las articulaciones es numérico
    if isnumeric(theta)

        % Crea una matriz NUMÉRICA vacía
        J = zeros(3, n);

    % De lo contrario, son simbólicos
    else

        % Crea una matriz SIMBÓLICA vacía
        J = sym(zeros(3, n));

    end

    % Itera a través de cada articulación
    for j = 1 : n

        % Si la articulación actual está en un marco de referencia que no afecta el marco analizado
        if j >= frame

            % Stop
            break;

        end

        % Cinemática directa hasta la j-ésima articulación
        H = forwardKinematicsHTM(DH, j + 1);

        % Actuation axis with respect to inertial frame
        z = H(1 : 3, 3);

        % Jacobian Matrix for the Angular Velocity
        J(:, j) = z;

    end
end