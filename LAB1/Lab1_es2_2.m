clear;
close all;
clc;

%% Es. 2.1 (Problema di Poisson a condizione miste con diff. centrate e ghost node)

fprintf("\n--- Es. 2.1 ---\n")

% Dati del problema
L = 1;
du0dx = 2;
uL = 1;
f = @(x) 4*pi^2*cos(2*pi*x);

% Discretizzazione scelta
N = 50;


% Calcolo della soluzione numerica
[x, uh] = poisson_misto_centrato(L, N, du0dx, uL, f);


% Soluzione esatta (per confronto)
uex = @(x) cos(2*pi*x) + 2*x - 2;


% Confronto grafico
figure
hold on

xplot = linspace(0, L, 1000);
plot(xplot, uex(xplot), 'linewidth', 2);

plot(x, uh, 'or', 'MarkerFaceColor', 'r');

title('Problema di Poisson (Neumann-Dirichlet b.c.) - NODO FANTASMA')
xlabel('x')
legend('u_{ex}', 'u_h', 'Location', 'northwest');

errore = max(abs(uh - uex(x)));
fprintf("Errore: %.2e.\n", errore)


%% Es. 2.2 (Rate di convergenza)

fprintf("\n--- Es. 2.2 ---\n")

NN = [50, 100, 200, 400, 800];

errori_max = zeros(length(NN), 1);
errori_h = zeros(length(NN), 1);

for i = 1:length(NN) 
    N = NN(i);
    h = L/N;

    [x, uh] = poisson_misto_centrato(L, N, du0dx, uL, f);

    errori_max(i) = max(abs(uh - uex(x)));
    errori_h(i) = sqrt(h)*norm(uh - uex(x));
end

p_max = log2(errori_max(1:end-1) ./ errori_max(2:end));
p_h = log2(errori_h(1:end-1) ./ errori_h(2:end));

fprintf("p_max:\n");
disp(p_max);

fprintf("p_h:\n");
disp(p_h);