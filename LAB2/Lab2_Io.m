%% Esercizio 1
clear;
clc;

% Dati iniziali
T = 3;
L = pi;
f = @(x,t) (2*sin(3.*x)-1).*exp(-t) + (x-pi)./(2*pi*(1+t).^2);
u0 = @(x) 2 + cos(x) +0.25*sin(3*x);
ua = @(t) 1;
ub = @(t) 1 + (t)./(1+t);

uex = @(x,t) (1 + cos(x) +0.25*sin(3*x)).*exp(-t) + (t.*(x+pi))./(2*pi.*(1+t)) + 1;

% 1 già risolto in altri file

% 2
N = 35;
K = 200;

[x,t,U] = calore_EA(L,N,T,K,ua,ub,f,u0);

% Supponendo che [x,t,U] = calore_EA(...) 
% ti restituisca
%   x: vettore colonna di dimensione (N+1)×1
%   t: vettore colonna di dimensione (K+1)×1
%   U: matrice (N+1)×(K+1) con U(i,j)=u(x_i,t_j)

% 1) Costruisci le griglie X e T
[X,T] = meshgrid(x,t);  
% meshgrid con (x,t) dà matrici di dimensione (K+1)×(N+1)

% 2) Plot della soluzione numerica EA
figure
surf(X, T, U.')   % trasponi U per fare corrispondere (row= t, col= x)
shading interp
title('Soluzione numerica (EA)')
xlabel('x')
ylabel('t')
zlabel('U(x,t)')
view(45,30)       % regola l’angolo di vista a piacere

% 3) Calcola la soluzione esatta sui nodi
Uex_mat = uex(X, T);   % stessa dimensione di X e T

% 4) Plot della soluzione esatta
figure
surf(X, T, Uex_mat)
shading interp
title('Soluzione esatta')
xlabel('x')
ylabel('t')
zlabel('u_{ex}(x,t)')
view(45,30)

% 5) (Opzionale) Confronto affiancato
figure
subplot(1,2,1)
surf(X, T, U.')
shading interp
title('EA')
xlabel('x'), ylabel('t'), zlabel('U')
view(45,30)

subplot(1,2,2)
surf(X, T, Uex_mat)
shading interp
title('Esatta')
xlabel('x'), ylabel('t'), zlabel('u_{ex}')
view(45,30)
