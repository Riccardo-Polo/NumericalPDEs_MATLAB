% MAIN PROGRAM
% risoluzione numerica dell'equazione del calore
% con condizioni di Dirichlet al bordo 

close all
clc

% Dati del problema
L = pi;
T = 3;
ua = @(t) 1; 
ub = @(t) (t./(1+t)) + 1; 
u0 = @(x) 2 + cos(x) + 0.25*sin(3*x); 
f  = @(x,t) (2*sin(3*x)-1)*exp(-t) + (x+pi)/(2*pi)/((1+t).^2);

% Soluzione esatta
uex = @(x,t) (1+cos(x)+0.25*sin(3*x)).*exp(-t) + 1 + (t./(1+t)).*(x+pi)/(2*pi);

% Discretizzazione e soluzione numerica
N = 35;
K = 200;
[x,t,u] = calore_CN(L,N,T,K,ua,ub,f,u0);


%% >> Visualizzazioni varie

% --- Animazione (plot dinamico) ---
for k = 1:K+1
    % Soluzione al tempo t_k
    plot(x, u(:, k));
    
    % Titolo interattivo + fix grafici
    title(sprintf('t = %.2f', t(k))); 
    xlabel('x');
    ylabel('u');
    axis([-L, L, 0, 3]); 
    
    pause(1.0/24.0);
end

% --- Grafico 3D (plot statico sul dominio spazio-tempo) ---

% Soluzione numerica
ax1 = subplot(1,2,1);
mesh(t, x, u)
title("Soluzione numerica")

% Soluzione esatta
ax2 = subplot(1,2,2);
space = linspace(-L, L, 200);
time = linspace(0, T, 100);
[xx, tt] = ndgrid(space, time);
mesh(tt, xx, uex(xx, tt));
title("Soluzione esatta")

hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});
rotate3d on