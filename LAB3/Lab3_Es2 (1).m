% Es. 2. Risoluzione numerica del problema di trasporto 1D con condizione 
% iniziale e condizione di inflow. Metodo a scelta tra: EA/C, UW ed EI/C.

clear;
close all;
clc;

% Dati del problema
a = -3; 
b = 3; 
T = 2;
c = 1;

u0 = %todo
f =  %todo
g =  %todo

% Soluzione esatta
uex = %todo


% Calcolo della soluzione numerica
N = 60;
K = %todo

[x, t, u_EA] = conservazione_EAC(a, b, N, T, K, c, f, g, u0);
[x, t, u_UW] = conservazione_UW( a, b, N, T, K, c, f, g, u0);
[x, t, u_EI] = conservazione_EIC(a, b, N, T, K, c, f, g, u0);

% Calcolo CFL
%...
CFL = %todo


% Rappresentazione della soluzione numerica confrontata con quella esatta
xesatta = linspace(a, b, 1000);

figure('units','normalized','outerposition',[0 0.25 1 0.5])
for k=1:K+1
 subplot(1, 4, 1)
 plot(xesatta, uex(xesatta, t(k)))
 axis([a, b, -1, 1])
 title(sprintf('Soluzione esatta [t = %.2f]', t(k)))
 xlabel('x')
 ylabel('u(x,t)')
 
 subplot(1, 4, 2)
 plot(x, u_EA(:, k), '-or')
 axis([a, b, -1, 1])
 title(sprintf('EA/C (CFL = %.2f)', CFL))
 xlabel('x')
 ylabel('u_h(x,t)')

 subplot(1, 4, 3)
 plot(x, u_UW(:, k), '-o', 'Color', [0.5, 0, 0.5])
 axis([a, b, -1, 1])
 title(sprintf('UW (CFL = %.2f)', CFL))
 xlabel('x')
 ylabel('u_h(x,t)')

 subplot(1, 4, 4)
 plot(x, u_EI(:, k), '-ok')
 axis([a, b, -1, 1])
 title(sprintf('EI/C (CFL = %.2f)', CFL))
 xlabel('x')
 ylabel('u_h(x,t)')
 
 pause(2.0/24.0);
end
