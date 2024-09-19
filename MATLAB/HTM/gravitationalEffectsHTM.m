function G = gravitationalEffectsHTM(m, g, DH, theta)
    
%% Esta función calcula los efectos gravitacionales del robot utilizando matrices de transformación homogéneas

%{
m: matriz con los valores de la masa de cada cuerpo rígido
g: vector de aceleración gravitacional (por ejemplo, [0 0 -9.80665]^T)
DH: parámetros de Denavit Hartenberg como una MATRIZ (simbólica o numérica)
theta: theta
%}

    % Obtener el número de cuerpos rígidos
    [rb, ~] = size(m);
    
    % Obtener el número de articulaciones
    [n, ~] = size(theta);

    % Si el vector de las articulaciones es numérico
    if isnumeric(theta)

        % Crea una matriz NUMÉRICA vacía.
        G = zeros(n, 1);

    % De lo contrario, son simbólicos
    else

        % Crea una matriz SIMBÓLICA vacía
        G = sym(zeros(n, 1));

    end

    % Iterar a través de todos los cuerpos rígidos
    for j = 1 : rb

        % Calcular la matriz jacobiana del cuerpo rígido actual
        J = jacobianMatrixHTM(DH, theta, j + 1);

        % Calcular la derivada actual de la energía potencia
        G = G + (m(j) * transpose(J) * g);

    end

end