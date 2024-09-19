function R = Rz(theta)

% Matriz de rotación con respecto al eje «z»
    R = [cos(theta) -sin(theta) 0 0
         sin(theta) +cos(theta) 0 0
         0.00000000 0.000000000 1 0
         0.00000000 0.000000000 0 1];
end