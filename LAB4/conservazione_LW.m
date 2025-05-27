function [x, t, U] = conservazione_LW(a, b, N, T, K, c, f, dfdt, dfdx, g, u0)
% ----  Risoluzione dell'equazione di conservazione scalare ----
%   u_t + c u_x = f nell'intervallo [a, b] 
%   con condizioni di inflow nel nodo sinistro e condizione iniziale.
%   Metodo numerico: LAX-WENDROFF.
% -----------------------------------------------
% Sintassi:
%   [x, t, U] = conservazione_LW(a, b, N, T, K, c, f, g, u0)
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

% Ciclo temporale
for k=1:K

    % Avanzamento nei nodi interni
    U(2:end-1, k+1) = U(2:end-1, k) + tau*f(x(2:end-1), t(k)) ... 
                                    - 0.5*lambda*(U(3:end, k) - U(1:end-2, k)) ...
                                    + 0.5*(lambda^2)*(U(3:end, k) - 2*U(2:end-1, k) + U(1:end-2, k)) ...
                                    + 0.5*(tau^2)*(dfdt(x(2:end-1), t(k)) - c*dfdx(x(2:end-1), t(k)));

    % Avanzamento nell'ultimo nodo
    U(end, k+1)     = U(end, k)     + tau*f(x(end), t(k)) ... 
                                    - lambda*(U(end, k) - U(end-1, k));
end
