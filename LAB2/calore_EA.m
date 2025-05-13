function [x, t, U] = calore_EA(L, N, T, K, ua, ub, f, u0)
% ----  Risoluzione dell'equazione del calore ----
%   u_t - u_xx = f nell'intervallo [-L,L] 
%   con condizioni al bordo di Dirichlet
%   e condizioni iniziali
%
%   Metodo numerico in spazio: DIFFERENZE FINITE CENTRATE.
%   Metodo numerico in tempo:  EULERO IN AVANTI.
% -----------------------------------------------
% Sintassi:
%   [x, t, U] = calore_EA(L, N, T, K, ua, ub, f, u0)
%
% Input:
%   L  (float)           semiampiezza intervallo spaziale (-L,L) 
%	N  (int)             numero di sottointervalli in (-L,L)
%   T  (float)           estremo finale intervallo temporale (0,T)
%	K  (int)             numero di sottointervalli in (0,T)
%   ua (handle function) condizione di Dirichlet in x=-L (in funzione di t)
%   ub (handle function) condizione di Dirichlet in x=L  (in funzione di t)
%   f  (handle function) funzione descrivente il termine noto
%   u0 (handle function) funzione descrivente la condizione iniziale (t=0)
%
% Output:
%   x (vettore N+1)       griglia spaziale
%   t (vettore K+1)       nodi temporali
%   U (matrice N+1 x K+1) soluzione numerica del problema

% Discretizzazione in spazio
h   = 2*L/N;
x = linspace(-L, L, N+1);
x = x';

% Discretizzazione in tempo
tau = T/K; 
t = linspace(0, T, K+1);
t = t';

% Inizializzazione della matrice soluzione U
U = zeros(N+1, K+1);

% Condizione iniziale
U(:, 1) = u0(x);

% Condizioni al bordo
U(1, :)   = ua(t);
U(end, :) = ub(t);

% Costruzione in formato sparso della matrice A
A = spdiags([-1 2 -1], [-1 0 1], N-1, N-1);
A = (1/h^2)*A;

% Matrice identità sparsa
I = speye(N-1, N-1);

% Ciclo iterativo
for k = 1:K
    
    % Assemblaggio termine noto al t_{k}
    F = f(x(2:end-1), t(k));
    
    % Correzione del termine noto con le condizioni al bordo
    F(1)  = F(1)   + ua(t(k))*(1/h^2);
    F(end)= F(end) + ub(t(k))*(1/h^2);

    % Calcolo della soluzione al tempo t_{k+1}
    U(2:end-1, k+1) = tau*F + (I - tau*A)*U(2:end-1, k);
end
