function R = Ry(beta)

% Matriz de rotación con respecto al eje «y»
    R = [+cos(beta) 0 sin(beta) 0
         0.00000000 1 0.000000 0
         -sin(beta) 0 cos(beta) 0
         0.00000000 0 0.000000 1];
end