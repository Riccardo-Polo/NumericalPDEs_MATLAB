% Analisi degli autovalori della matrice A=1/h^2 * tridiag(-1,2,-1)
% associata alla discretizzazione in spazio di u_xx con differenze 
% finite del secondo ordine

clear;
close all;
clc;

disp('Es 1.4')
disp('Analisi degli autovalori della matrice A = 1/h^2 * tridiag(-1,2,-1)')
disp('associata alla discretizzazione in spazio di u_xx con differenze')
disp('finite del secondo ordine')

L = pi;
N = 100;
M = 5;

h = zeros(M, 1);
lambdamin = zeros(M, 1);
lambdamax = zeros(M, 1);

for i = 1:M
    % Costruzione in formato sparso della matrice A
    h(i) = 2*L/N;
    A = (1/h(i)^2) * spdiags([-1 2 -1], [-1 0 1], N-1, N-1);
    N = 2*N;
    % Calcolo degli autovalori minimo e massimo
    lambdas = eig(A);
    lambdamin(i) = min(lambdas);
    lambdamax(i) = max(lambdas);
end

loglog(h, lambdamin);
hold on 
loglog(h, lambdamax);
loglog(h, h.^-2, '--k');
xlabel('h')
legend('\lambda_{min}','\lambda_{max}','h^{-2}')