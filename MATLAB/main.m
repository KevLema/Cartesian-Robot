clc; clear variables

%% Cargar directorios
addpath("HTM/");
addpath("../ROS/");

%% Definir los atributos del robot (simbólicos o numéricos)

% Variables de las articulaciones
theta = sym("theta_", [3, 1]);
dtheta   = sym("dtheta_", [3, 1]);
ddtheta   = sym("ddtheta_", [3, 1]);

% Valores deseados de las articulaciones
dtheta_d = sym("dtheta_d_", [3, 1]);
ddtheta_d = sym("ddtheta_d_", [3, 1]);

% Variables del efector final
X     = [sym("x"); sym("y"); sym("z")];
dX   = [sym("dx");   sym("dy"); sym("dz")];
ddX   = [sym("ddx"); sym("ddy"); sym("ddz")];

% Valores deseados del efector final
dX_d = [sym("dx_d"); sym("dy_d"); sym("dz_d")];
ddX_d = [sym("ddx_d"); sym("ddy_d"); sym("ddz_d")];

% Parámetros físicos
g = [0 0 -9.80665]';
m = sym("m_", [3, 1]);
m_p = sym("m_p");
L = sym("l_", [4, 1]);
dt = 0.01;

% Parámetros de control
syms e_theta_1(t) e_theta_2(t) e_theta_3(t) e_x(t) e_y(t) e_z(t) t
e  = [ e_theta_1(t);  e_theta_2(t);  e_theta_3(t);  diff([e_theta_1(t);  e_theta_2(t);  e_theta_3(t)])];
de = diff(e);

e_c  = [e_x(t); e_y(t); e_z(t); diff([e_x(t); e_y(t); e_z(t)])];
de_c = diff(e_c);

K   = sym("k_", [6, 6]);
K_c = sym("k_c_", [6, 6]);

%% Cinemática del robot

% Parámetros de Denavit-Hartenberg
DH = denavitHartenberg(theta, L);

% Cinemática directa usando matrices de transformación homogéneas
H = forwardKinematicsHTM(DH, 5);

%% Cinemática diferencial del robot

% Matriz Jacobiana geométrica usando matrices de transformación homogéneas
J_g = jacobianMatrixHTM(DH, theta, 5);

% Verificar si se puede realizar la cinemática inversa numérica
if isnumeric(theta)
    
    % Cinemática inversa usando la matriz Jacobiana geométrica
    theta_htm = inverseKinematicsHTM(zeros(3, 1), L, H(1 : 3, 4), eye(3), 5, dt);
    
    % Cinemática inversa usando matrices de transformación homogéneas
    [n, ~] = size(theta);
    legends = string(zeros(1, n));
    figure()
    hold on;
    for i = 1 : n
        plot(theta_htm(i, :), 'LineWidth', 1)
        legends(1, i) = string(append("$\theta_", num2str(i), "$"));
    end
    title("Inverse Kinematics using Geometric Jacobian Matrix", 'Interpreter', 'latex')
    xlabel("Samples $\left[ k \right]$", 'Interpreter', 'latex')
    ylabel("Amplitude $\left[ meters \right]$", 'Interpreter', 'latex')
    legend(legends, 'Interpreter', 'latex', 'FontSize', 12)
    grid on;
    hold off;

end

%% Dinámica del robot

% Invertible por la izquierda:
% filas > columnas
% 
% Invertible por la derecha:
% filas < columnas
% 
% Invertible por ambos lados:
% filas = columnas

% Matriz de masa usando matrices de transformación homogéneas
D = massMatrixHTM(m, DH, theta);

% Matriz de efectos centrífugos y de Coriolis usando matrices de transformación homogéneas
C = centrifugalCoriolis(D, theta, dtheta);

% Efectos gravitacionales usando matrices de transformación homogéneas
G = gravitationalEffectsHTM(m, g, DH, theta);

% Torques generados por las placas
tau_p = m_p * transpose(J_g) * g;

% Modelo dinámico en el espacio de estados
A = [zeros(3, 3) eye(3)
     zeros(3, 3) -D \ C];
x = [theta; dtheta];
Phi = [zeros(3, 1); D \ (tau_p - G)];
B = [zeros(3, 3) zeros(3, 3)
     zeros(3, 3) inv(D)];
u = [zeros(3, 1); sym("tau_", [3, 1])];

% Model
dx = (A * x) + Phi + (B * u);

%% Dinámica del error de las articulaciones

% Control de articulaciones
u = pinv(B) * ([dtheta_d; ddtheta_d] + (K * e) - (A * x) - Phi);

% Dinámica del error
f_e = subs(de - [dtheta_d; ddtheta_d] + (A * x) + Phi + (B * u), [dtheta(1) - dtheta_d(1); dtheta(2) - dtheta_d(2); dtheta(3) - dtheta_d(3)], -diff(e(1 : 3, :)));

% Simplificación de la dinámica del error
f_e = subs(f_e, [K(4, 2), K(4, 3), K(4, 5), K(4, 6), K(5, 1), K(5, 3), K(5, 4), K(5, 6), K(6, 1), K(6, 2), K(6, 4), K(6, 5)], zeros(1, 12)) == 0;

%% Control

% Parámetros físicos
m_p = 20;
m   = [7.538; 4.211; 1.202];
L   = [0.755 0.1269 0.083 0.083]';

% Parámetros de las articulaciones
theta_0  = zeros(3, 1);
dtheta_0 = zeros(3, 1);

%% Ganancias

% % Subamortiguado
% k_d = 2;
% k_p = (k_d^2 / 4) + 7;

% Críticamente amortiguado
k_d = 7;
k_p = (k_d^2 / 4);

% % Sobreamortiguado
% k_d = 2;
% k_p = (k_d^2 / 4) - 1;

% Matrices de ganancias
K   = [zeros(3, 3) zeros(3, 3)
         k_p  0  0     k_d 0  0
           0 k_p 0      0 k_d 0
           0  0 k_p     0  0 k_d];

% K_c = [zeros(3, 3) zeros(3, 3)
%          k_p  0  0     k_d 0  0
%            0 k_p 0      0 k_d 0
%            0  0 k_p     0  0 k_d];

%% Planificación de trayectorias (para las articulaciones)

% Parámetros para las placas
d   = 2.54 / 200;                       % Plate's thickness. 1/2" to meters
d_r = 1 / 100;                          % Base de la placa. 1 centímetro a metros
h_s = 1 / 100;                          % Altura de seguridad
n   = 7;                                % Número de placas

% Coordenadas euclidianas
x_1 = 0.2;
y_1 = 0.2;
z_1 = 0.2;

% Coordenadas euclidianas
x_2 = 0.4;
y_2 = 0.4;

% Inicializar variable para almacenar las trayectorias
theta_d = theta_0;
theta_d_gz = zeros(3, 7 * n);

% Itera a través de todas las placas
for i = 1 : n

    % Punto de partida: (x_1, y_1, z_1)
    theta_d                  = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_1 y_1 z_1]', eye(3), 5, 0.001)];    % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 6) = theta_d(:, end);                                                                               % For Gazebo (partial trajectory, only keypoints)

    % Bajar para una placa: (x_1, y_1, z_1 - i * (d + d_r) - h_s)
    theta_d = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_1 y_1 z_1 - i * (d + d_r) - h_s]', eye(3), 5, 0.001)]; % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 5) = theta_d(:, end);                                                                                 % For Gazebo (partial trajectory, only keypoints)
    
    %{
        *** EL AGENTE DEBE SER ACTIVADO EN ESTE PUNTO PARA TRANSPORTAR LA i-ésima PLACAE ***
    %}

    % Subir con la placa: (x_1, y_1, z_1)
    theta_d = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_1 y_1 z_1]', eye(3), 5, 0.001)];   % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 4) = theta_d(:, end);                                                             % For Gazebo (partial trajectory, only keypoints)

    % Moverse a la parte superior de la columna de placas: (x_2, y_2, z_1)
    theta_d = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_2 y_2 z_1]', eye(3), 5, 0.001)];   % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 3) = theta_d(:, end);                                                             % For Gazebo (partial trajectory, only keypoints)

    % Bajar para alcanzar el punto más cercano de descarga de las placas: (x_2, y_2, z_1 + (i - n) * (d + d_r) + h_s)
    theta_d = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_2 y_2 z_1 + (i - n) * (d + d_r) + h_s]', eye(3), 5, 0.001)];   % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 2) = theta_d(:, end);                                                                                         % For Gazebo (partial trajectory, only keypoints)

    % Alcanzar el punto de descarga de las placas: (x_2, y_2, z_1 + (i - n) * (d + d_r))
    theta_d = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_2 y_2 z_1 + (i - n) * (d + d_r)]', eye(3), 5, 0.001)];         % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 1) = theta_d(:, end);                                                                                         % For Gazebo (partial trajectory, only keypoints)

    %{
        *** EL IMÁN DEBE SER DESACTIVADO EN ESTE PUNTO PARA DESCARGAR LA i-ésima PLACA ***
    %}

    % Moverse a la parte superior de la columna de placas: (x_2, y_2, z_1)
    theta_d = [theta_d(:, :) inverseKinematicsHTM(theta_d(:, end), L, [x_2 y_2 z_1]', eye(3), 5, 0.001)];   % For the controller (full trajectory)
    theta_d_gz(:, 7 * i - 0) = theta_d(:, end);                                                             % For Gazebo (partial trajectory, only keypoints)

end

% Calcular las derivadas de las trayectorias
dtheta_d  = gradient(theta_d);
ddtheta_d = gradient(dtheta_d);

% Guardar el comando para Gazebo con los puntos calculados
file = fopen('../ROS/gazebo_command.txt', 'w');
fprintf(file, append('gz topic -t "/model/euclidian_robot/joint_trajectory" -m gz.msgs.JointTrajectory -p ', "'joint_names: ", '"joint_1"; joint_names: "joint_2"; joint_names: "joint_3";'));

% Iterar a través de todos los puntos
for i = 1 : 7 * n

    % Agregar datos al comando
    fprintf(file, "points { positions: %f; positions: %f; positions: %f; time_from_start {  sec: %d; nsec: 0 } }; ",  theta_d_gz(1, i), theta_d_gz(2, i), theta_d_gz(3, i) - 1.1075, i);

end

fprintf(file, "'");

% Cerrar archivo
fclose(file);

%% Preliminares de la simulación

% Tiempo de simulación (segundos)
t = 1500;

% Número de muestras para la simulación
[~, k] = size(theta_d);

% Transformar en series temporales
theta_d   = timeseries( theta_d, linspace(0, t, k));
dtheta_d  = timeseries(dtheta_d, linspace(0, t, k));
ddtheta_d = timeseries(ddtheta_d, linspace(0, t, k));

%% Simulación (modelo de articulaciones)

% Simulate
sim("joints_2022b.slx");

%% Gráficas

% Posiciones de las articulaciones
figure()
hold on;
plot(x.time, reshape(x.data(1, :, :), 1, []), 'LineWidth', 1)
plot(x.time, reshape(x.data(2, :, :), 1, []), 'LineWidth', 1)
plot(x.time, reshape(x.data(3, :, :), 1, []), 'LineWidth', 1)
plot(theta_d, '--', 'LineWidth', 1)
title("Position of each Joint", 'Interpreter', 'latex')
xlabel("Time $\left[ sec \right]$", 'Interpreter', 'latex')
ylabel("Amplitude $\left[ meters \right]$", 'Interpreter', 'latex')
legend(["$\theta_1 \left( t \right)$", "$\theta_2 \left( t \right)$", "$\theta_3 \left( t \right)$", "$\theta_{d_1} \left( t \right)$", "$\theta_{d_2} \left( t \right)$", "$\theta_{d_3} \left( t \right)$"], 'Interpreter', 'latex', 'FontSize', 12)
grid on;
hold off;

% Velocidades de las articulaciones
figure()
hold on;
plot(x.time, reshape(x.data(4, :, :), 1, []), 'LineWidth', 1)
plot(x.time, reshape(x.data(5, :, :), 1, []), 'LineWidth', 1)
plot(x.time, reshape(x.data(6, :, :), 1, []), 'LineWidth', 1)
plot(dtheta_d, '--', 'LineWidth', 1)
title("Velocity of each Joint", 'Interpreter', 'latex')
xlabel("Time $\left[ sec \right]$", 'Interpreter', 'latex')
ylabel("Amplitude $\left[ \frac{m}{s} \right]$", 'Interpreter', 'latex')
legend(["$\dot{\theta}_1 \left( t \right)$", "$\dot{\theta}_2 \left( t \right)$", "$\dot{\theta}_3 \left( t \right)$", "$\dot{\theta}_{d_1} \left( t \right)$", "$\dot{\theta}_{d_2} \left( t \right)$", "$\dot{\theta}_{d_3} \left( t \right)$"], 'Interpreter', 'latex', 'FontSize', 12)
grid on;
hold off;

% Fuerzas del efector final
figure()
hold on;
plot(force_c.time, reshape(force_c.data(1, :, :), 1, []), '-r', 'LineWidth', 1)
plot(force_c.time, reshape(force_c.data(2, :, :), 1, []), '-g', 'LineWidth', 1)
plot(force_c.time, reshape(force_c.data(3, :, :), 1, []), '-b', 'LineWidth', 1)
title("Force at the End Effector", 'Interpreter', 'latex')
xlabel("Time $\left[ sec \right]$", 'Interpreter', 'latex')
ylabel("Amplitude $\left[ Newtons \right]$", 'Interpreter', 'latex')
legend(["$f_x \left( t \right)$", "$f_y \left( t \right)$", "$f_z \left( t \right)$"], 'Interpreter', 'latex', 'FontSize', 12)
grid on;
hold off;

% Función de error
figure()
plot(e, 'LineWidth', 1)
title("Error Function", 'Interpreter', 'latex')
xlabel("Time $\left[ sec \right]$", 'Interpreter', 'latex')
ylabel("Amplitude", 'Interpreter', 'latex')
legend(["$e_{\theta_1} \left( t \right)$", "$e_{\theta_2} \left( t \right)$", "$e_{\theta_3} \left( t \right)$", "$\dot{e}_{\theta_1} \left( t \right)$", "$\dot{e}_{\theta_2} \left( t \right)$", "$\dot{e}_{\theta_3} \left( t \right)$"], 'Interpreter', 'latex', 'FontSize', 12)
grid on;

% Función de control
figure()
plot(force, 'LineWidth', 1)
title("Control Function", 'Interpreter', 'latex')
xlabel("Time $\left[ sec \right]$", 'Interpreter', 'latex')
ylabel("Amplitude $\left[ Newtons \right]$", 'Interpreter', 'latex')
legend(["$f_{\theta_1} \left( t \right)$", "$f_{\theta_2} \left( t \right)$", "$f_{\theta_3} \left( t \right)$"], 'Interpreter', 'latex', 'FontSize', 12)
grid on;