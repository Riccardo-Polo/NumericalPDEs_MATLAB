% Valutazione dell'ordine di convergenza del metodo UPWIND
% applicato al problema di trasporto (leggi di conservazione)

clear;
close all;
clc;

disp('Es 2')
disp('Valutazione dell ordine di convergenza (metodo a scelta tra LF ed LW)')
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
e1s = zeros(M, 1);
e2s = zeros(M, 1);

for i = 1:M
    %[x, t, uh] = conservazione_LF(a, b, N, T, K, c, f, g, u0);
    [x, t, uh] = conservazione_LW(a, b, N, T, K, c, f, dfdt, dfdx, g, u0);
    h = (b-a)/N;

    [xx, tt] = ndgrid(x, t);
    e1s(i) = max(max(abs(uh-uex(xx, tt))));
    e2s(i) = max(sqrt(h*sum((uh-uex(xx, tt)).^2)));

    N = 2*N;
    K = 2*K;
end

p1 = log2(e1s(1:end-1)./e1s(2:end));
p2 = log2(e2s(1:end-1)./e2s(2:end));

p1
p2