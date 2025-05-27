a = 0.022; % dispersion prefactor
N = 256;
tf = 2;
dt = 0.001;
Nt = round(tf/dt);
x0 = -pi;
xf = pi;
x = (xf-x0)/N*(-N/2:N/2-1);
u = cos(x); % initial condition
uu = zeros(Nt,N); uu(1,:) = u;

M = round(N/2);
k = zeros(1,N);
k(1:M) = 0:M-1; k(M+1:N) = -M:-1;
ik3 = 1j*k.^3;
U = fft(u); % go to spectral space
for n = 2:Nt
    t = n*dt;
    g = -0.5j*dt*k;
    expaik = exp(a*dt*ik3/2); % integrating
    expaik2 = expaik.^2; % factor
    % Runge-Kutta scheme
    rk1 = g.*fft(real(ifft(U)).^2);
    rk2 = g.*fft(real(ifft(expaik.*(U+rk1/2))).^2);
    rk3 = g.*fft(real(ifft(expaik.*U+rk2/2)).^2);
    rk4 = g.*fft(real(ifft(expaik2.*U+expaik.*rk3)).^2);
    U = expaik2.*U + (expaik2.*rk1 + 2*expaik.*(rk2+rk3) + rk4)/6;
    u = real(ifft(U)); % back to real space
    uu(n,:) = u;
end


figure
for k=1:10:Nt
    plot(x', uu(k,:), x', uu(1,:), "-r" )
    axis([-3,3,-1,3])
    pause(0.05);
    
end


