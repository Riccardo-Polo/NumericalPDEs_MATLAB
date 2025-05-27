% Valutazione dell'ordine di convergenza del metodo UPWIND
% applicato al problema di trasporto (leggi di conservazione)

clear;
close all;
clc;

disp('Es 2')
disp('Valutazione dell ordine di convergenza (LF vs LW)')
disp('applicato al problema di trasporto (leggi di conservazione)')

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

% Stima dell'ordine per dimezzamenti successivi
N = 300;
K = 200;
M = 5;

% TODO
% TODO
% TODO
% TODO