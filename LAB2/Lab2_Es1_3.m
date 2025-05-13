% Valutazione dell'ordine di convergenza relativi a schemi EA, EI e CN per
% l'equazione del calore

clear;
close all;
clc;

disp('Es 1.3')
disp('Valutazione dell''ordine di convergenza per l''equazione del calore')

% Dati del problema
L = pi;
T = 3;
ua = @(t) 1; 
ub = @(t) (t./(1+t)) + 1; 
u0 = @(x) 2 + cos(x) + 0.25*sin(3*x); 
f  = @(x,t) (2*sin(3*x)-1)*exp(-t) + (x+pi)/(2*pi)/((1+t).^2);

% Soluzione esatta
uex = @(x,t) (1+cos(x)+0.25*sin(3*x)).*exp(-t) + 1 + (t./(1+t)).*(x+pi)/(2*pi);

% Valutazione dell'ordine di convergenza p (risp. errori e1 ed e2)
M = 4;
e1s = zeros(M, 1);
e2s = zeros(M, 1);

% Ciclo sui dimezzamenti
N = 25; 
K = 35;
for i = 1:M
    [x, t, u] = calore_EI(L, N, T, K, ua, ub, f, u0);
    h = 2*L/N;

    [xx,tt] = ndgrid(x,t);
    uEX = uex(xx,tt);

    % Calcolo errori
    e1s(i) = max(max(abs(u-uEX)));
    e2s(i) = max((sqrt(h*sum((u-uEX).^2))));
    
    % Dimezzamento
    K = 2*K;
    N = 2*N;
end


% Stima dell'ordine di convergenza
p1 = log2(e1s(1:end-1)./e1s(2:end)) ;
p2 = log2(e2s(1:end-1)./e2s(2:end)) ;

fprintf("\np1 =\n")
disp(p1)

fprintf("\np2 =\n")
disp(p2)