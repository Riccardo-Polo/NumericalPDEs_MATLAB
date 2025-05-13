function [x, u] = poisson_misto_centrato(L, N, du0dx, uL, f)
% ----  Risoluzione dell'equazione modello ----
%   -u''= f   con condizioni al bordo miste Neumann-Dirichlet 
%             u'(0)=du0dx, u(L)=uL nell'intervallo [0,L]
%
% Metodo numerico: differenze finite centrate.
% -----------------------------------------------
% Sintassi:
%   [x, u] = poisson_misto_centrato(L, N, du0dx, uL, f)
%
% Input:
%   L     (float)           lunghezza dell'intervallo [0, L]
%   N     (int)             numero di sottointervalli in (0, L)
%   du0dx (float)           condizione di Neumann in x = 0
%   uL    (float)           condizione di Dirichlet in x = L
%   f     (function handle) funzione descrivente la forzante / termine noto
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
e = ones(5*N, 1);
A = spdiags([-e 2*e -e],[-1 0 1], N, N);

% Correzione della matrice A (ghost node vicino a x = 0)
A(1, 1) = 1;
A(1, 2) = -1;

% Costruzione del termine noto F
F = f(x(1:end-1));
F = F*(h^2);

% Correzione del termine noto (inclusione delle condizioni al bordo)
F(1)   = 0.5*F(1)   - h*du0dx; % Neumann (ghost)
F(end) = F(end) + uL;          % Dirichlet

% Risoluzione del sistema lineare
u = A\F;
u = [u; uL];
 

