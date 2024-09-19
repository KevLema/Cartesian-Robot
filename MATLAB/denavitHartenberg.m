function DH = denavitHartenberg(theta, L)

%% Cartesian Robot
    DH = [0 0 0 0
          0.000 L(1)            0 -pi/2
          -pi/2 L(2) + theta(1) 0 -pi/2
          +pi/2 L(3) + theta(2) 0 -pi/2
              0 L(4) + theta(3) 0 0.00];
    
end
