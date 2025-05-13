clear;
clc;

%% ESERCIZIO 1
clear;
clc;

fprintf("\n--- Es. 1.1 ---\n")
fprintf("\n Vedi file poisson_dirichlet_centrato\n")

fprintf("\n--- Es. 1.2 ---\n")

% Dichiaro i dati forniti
L = 1;
f = @(x) 96*x.^2;
u0 = 0;
uL = -1;
N = 40;
uex = @(x) 7*x - 8*x.^4;

% Calcolo soluzione numerica
[x,u] = poisson_dirichlet_centrato(L,N,u0,uL,f);


% Confronto grafico
figure
hold on

xplot = linspace(0, L, 1000);
plot(xplot, uex(xplot));

plot(x, u, 'or', 'MarkerFaceColor', 'r');

title('Problema di Poisson (Dirichlet b.c.)')
xlabel('x')
legend('u_{ex}', 'u')

errore = max(abs(u - uex(x)));
fprintf("Errore: %.2e.\n", errore)

fprintf("\n--- Es. 1.3 ---\n")
Ni = [50,100,200,400,800];
errinf = ones(length(Ni),1);
errh = ones(length(Ni),1);

pmax = ones(length(Ni),1);
ph = ones(length(Ni),1);

for i = 1:length(Ni)
    N = Ni(i);
    h = L/N;
    [x,u] = poisson_dirichlet_centrato(L,N,u0,uL,f);

    errinf(i,1) = max(abs(u-uex(x)));
    errh(i,1) = sqrt(h)*norm(u-uex(x));
end

p_max = log2(errinf(1:end-1) ./ errinf(2:end));
p_h = log2(errh(1:end-1) ./ errh(2:end));

fprintf("p_max:\n");
disp(p_max);

fprintf("p_h:\n");
disp(p_h);

