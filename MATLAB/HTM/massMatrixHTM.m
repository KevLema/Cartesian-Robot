function D = massMatrixHTM(m, DH, theta)
    
%% Esta función calcula la matriz de masa del robot 
%{
    m: array con los valores de la masa de cada cuerpo rígido
DH: parámetros de Denavit-Hartenberg como un ARREGLO (simbólico o numérico)
theta: theta
frame: el marco en movimiento, como un NÚMERO ENTERO, a estudiar (considerado como la n-ésima fila del arreglo mencionado arriba)
%}

    % Obtener el número de cuerpos rígidos
    [rb, ~] = size(m);
    
    % Obtener el número de articulaciones
    [n, ~] = size(theta);

    % Si el vector de articulaciones es numérico
    if isnumeric(theta)

        % Crea una matriz NUMÉRICA vacía
        D = zeros(n, n);

    % De lo contrario, son simbólicos
    else

        % Crea una matriz SIMBÓLICA vacía
        D = sym(zeros(n, n));

    end

    % Iterar a través de todos los cuerpos rígidos
    for j = 1 : rb
        
        % Calcular la matriz jacobiana del marco actual
        J = jacobianMatrixHTM(DH, theta, j + 1);

        % Calcular la matriz de inercia actual
        D = D + (m(j) * transpose(J) * J);

    end
    
end