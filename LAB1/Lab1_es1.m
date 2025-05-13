clear;
close all;
clc;

%% Es. 1.1 e 1.2 (Risoluzione del problema di Poisson per un dato N)

fprintf("\n--- Es. 1.1 e 1.2 ---\n")

% Dati del problema
L = 1;
u0 = 0;
uL = -1;
f = @(x) 96*x.^2;

% Discretizzazione scelta
N = 40;


% Calcolo della soluzione numerica
[x, uh] = poisson_dirichlet_centrato(L, N, u0, uL, f);


% Soluzione esatta (per confronto)
uex = @(x) -8*x.^4 +7*x;


% Confronto grafico
figure
hold on

xplot = linspace(0, L, 1000);
plot(xplot, uex(xplot), 'linewidth', 2);

plot(x, uh, 'or', 'MarkerFaceColor', 'r');

title('Problema di Poisson (Dirichlet b.c.)')
xlabel('x')
legend('u_{ex}', 'u_h')

errore = max(abs(uh - uex(x)));
fprintf("Errore: %.2e.\n", errore)


%% Es. 1.3 (Rate di convergenza)

fprintf("\n--- Es. 1.3 ---\n")

NN = [50, 100, 200, 400, 800];

errori_max = zeros(length(NN), 1);
errori_h = zeros(length(NN), 1);

for i = 1:length(NN) 
    N = NN(i);
    h = L/N;

    [x, uh] = poisson_dirichlet_centrato(L, N, u0, uL, f);

    errori_max(i) = max(abs(uh - uex(x)));
    errori_h(i) = sqrt(h)*norm(uh - uex(x));
end

p_max = log2(errori_max(1:end-1) ./ errori_max(2:end));
p_h = log2(errori_h(1:end-1) ./ errori_h(2:end));

fprintf("p_max:\n");
disp(p_max);

fprintf("p_h:\n");
disp(p_h);
