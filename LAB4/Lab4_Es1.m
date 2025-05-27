% Es. 2. Risoluzione numerica del problema di trasporto 1D con condizione 
% iniziale e condizione di inflow. Metodo a scelta tra: EA/C, UW ed EI/C.

clear;
close all;
clc;

disp('Es 1')
disp('Esempio di utilizzo della function per la risoluzione del problema')
disp('conservazione con condizioni di inflow e condizione iniziale (LF vs LW)')

% Dati del problema
a = -3; 
b = 3; 
T = 2;
c = 1;
u0 = @(x) cos(pi*x).^4.* (abs(x) <= 0.5);
f = @(x,t) 0*x.*t;
dfdt = @(x,t) 0*x.*t;
dfdx = @(x,t) 0*x.*t;
g = @(t) 0*t;

% Soluzione esatta
uex = @(x,t) u0(x-c*t);


% Calcolo della soluzione numerica
N = 60;
K = 40; % per CFL = 0.5; usare K = 10 per CFL = 2.

[x, t, u_LF] = conservazione_LF(a, b, N, T, K, c, f, g, u0);
[x, t, u_LW] = conservazione_LW(a, b, N, T, K, c, f, dfdt, dfdx, g, u0);

% Calcolo CFL
h = (b-a)/N;
tau = T/K;
CFL = c*tau/h;


% Rappresentazione della soluzione numerica confrontata con quella esatta
xesatta = linspace(a, b, 1000);

figure('units','normalized','outerposition',[0 0.25 1 0.5])
for k=1:K+1
 subplot(1, 3, 1)
 plot(xesatta, uex(xesatta, t(k)))
 axis([a, b, -1, 1])
 title(sprintf('Soluzione esatta [t = %.2f]', t(k)))
 xlabel('x')
 ylabel('u(x,t)')
 
 subplot(1, 3, 2)
 plot(x, u_LF(:, k), '-or')
 axis([a, b, -1, 1])
 title(sprintf('Lax-Friedrichs (CFL = %.2f)', CFL))
 xlabel('x')
 ylabel('u_h(x,t)')

 subplot(1, 3, 3)
 plot(x, u_LW(:, k), '-o', 'Color', [0.5, 0, 0.5])
 axis([a, b, -1, 1])
 title(sprintf('Lax-Wendroff (CFL = %.2f)', CFL))
 xlabel('x')
 ylabel('u_h(x,t)')
 
 pause(2.0/24.0);
end
