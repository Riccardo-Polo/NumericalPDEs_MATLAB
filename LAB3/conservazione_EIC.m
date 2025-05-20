function [x, t, U] = conservazione_EIC(a, b, N, T, K, c, f, g, u0)
% ----  Risoluzione dell'equazione di conservazione scalare ----
%   u_t + c u_x = f nell'intervallo [a, b] 
%   con condizioni di inflow nel nodo sinistro e condizione iniziale.
%   Metodo numerico: EI in tempo, DIFF. CENTRATE in spazio.
% -----------------------------------------------
% Sintassi:
%   [x, t, U] = conservazione_EIC(a, b, N, T, K, c, f, g, u0)
%
% Input:
%   a  (float)           estremo inferiore intervallo spaziale
%   b  (float)           estremo superiore intervallo spaziale
%	N  (int)             numero di sottointervalli in (a,b)
%   T  (float)           tempo finale
%	K  (int)             numero di sottointervalli in (0,T)
%   c  (float)           coefficiente di trasporto
%   f  (handle function) funzione descrivente il termine noto
%   g  (handle function) condizione di Dirichlet (inflow) in x = a,
%                        assegnata come funzione di t
%   u0 (handle function) funzione descrivente la condizione iniziale (t=0)
%
% Output:
%   x (vettore N+1)       griglia spaziale
%   t (vettore K+1)       nodi temporali
%   U (matrice N+1 x K+1) soluzione numerica del problema

% Discretizzazione in spazio e tempo
x = linspace(a, b, N+1)';
t = linspace(0, T, K+1)';
h   = (b-a)/N; 
tau = T/K;
lambda = c*tau/h;

% Inizializzazione della matrice soluzione u
U = zeros(N+1, K+1);

% Condizione iniziale
U(:, 1) = u0(x);  

% Inflow
U(1, 1:end) = g(t(1:end));

% Pre-assemblaggio matrice A
A = spdiags([-lambda/2, 1, lambda/2], [-1, 0, 1], N, N);
A(end, end) = 1+lambda;
A(end, end-1) = -lambda;

% Ciclo temporale
for k=1:K

    % Assemblaggio termine noto
    F = tau*f(x(2:end), t(k+1)) + U(2:end, k);
    F(1) = F(1) + 0.5*lambda*U(1, k+1);

    % Avanzamento temporale
    U(2:end, k+1) = A \ F ;
end
