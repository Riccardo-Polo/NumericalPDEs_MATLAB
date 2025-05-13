function [x, u] = poisson_dirichlet_centrato(L, N, u0, uL, f)
% ----  Risoluzione dell'equazione modello ----
%   -u''= f   con condizioni al bordo di Dirichlet u(0)=u0, u(L)=uL
%             nell'intervallo [0,L]
%
% Metodo numerico: differenze finite centrate.
% -----------------------------------------------
% Sintassi:
%   [x, u] = poisson_dirichlet_centrato(L, N, u0, uL, f)
%
% Input:
%   L  (float)           lunghezza dell'intervallo [0, L]
%   N  (int)             numero di sottointervalli in (0, L)
%   u0 (float)           condizione di Dirichlet in x = 0
%   uL (float)           condizione di Dirichlet in x = L
%   f  (function handle) funzione descrivente la forzante / termine noto
%
%
% Output:
%   x  (vettore col. N+1)   punti della griglia in cui viene approssimata
%                           la soluzione
%   u  (vettore col. N+1)   soluzione numerica del problema modello


% Costruzione della griglia
x = linspace(0, L, N+1)'; 
h = L/N;

% Costruzione della matrice A
e = ones(N-1, 1);
A = spdiags([-e 2*e -e],[-1 0 1], N-1, N-1);

% Costruzione del termine noto F
F = f(x(2:end-1));
F = F*(h^2);

% Correzione del termine noto (inclusione delle condizioni al bordo)
F(1)   = F(1)   + u0;
F(end) = F(end) + uL;

% Risoluzione del sistema lineare
u = A\F;
u = [u0; u; uL];
 

