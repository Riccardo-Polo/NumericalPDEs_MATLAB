function [x, t, U] = calore_EA(L, N, T, K, ua, ub, f, u0)
% ----  Risoluzione dell'equazione del calore ----
%   u_t - u_xx = f nell'intervallo [-L,L] 
%   con condizioni al bordo di Dirichlet
%   e condizioni iniziali
%
%   Metodo numerico in spazio: DIFFERENZE FINITE CENTRATE.
%   Metodo numerico in tempo:  XXX.
% -----------------------------------------------
% Sintassi:
%   [x, t, U] = calore_XXX(L, N, T, K, ua, ub, f, u0)
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

% costruzione griglia
x = linspace(-L,L,N+1);
h = 2*L/N;

% Costruzione matrice A
e = ones(N,1);
A = spdiags([-e 2*e -e], -1:1, N, N);
F = ones(N,1).*f(x);

% Discretizzazione in tempo
% ...
% ...

% Inizializzazione della matrice soluzione U
% ...
% ...

% Condizione iniziale
% ...
% ...

% Condizioni al bordo
% ...
% ...

% Costruzione in formato sparso della matrice A
% ...
% ...

% Matrice identità sparsa
I = speye(N-1, N-1);

% Ciclo iterativo
for k = 1:K
    
    % Assemblaggio termine noto
    % ...
    % ...
    
    % Correzione del termine noto con le condizioni al bordo
    % ...
    % ...

    % Calcolo soluzione al tempo t_{k+1}
    % ...
    % ...
    
end
