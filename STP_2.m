%% 1
clear
Q = 3;
b = 0.5;
T = 1;
M = 10^4;
N = 100;
TAU = 6; %v zadani indexovano od 0 (6 hodnot)
K = 95; %v zadani indexovano od 0 (95 hodnot)
sigma_bilysum = Q*(1-exp(-2*b*T));

% generovani realiazaci
Xk = zeros(N, M);
for i = 1:M
    Xk(1, i) = randn * sqrt(Q); % prvni X
    for j = 2:N
        % X_{k+1} = e^{-bT}  * X_k          + W_k
        Xk(j, i) = exp(-b*T) * Xk(j - 1, i) + randn * sqrt(sigma_bilysum);
    end
end

% odhad parametru
cov_Xk = zeros(TAU, K);
for k = 1:K
    for tau = 1:TAU    
        Xk_a = Xk(k, :);
        Xk_b = Xk(k + tau - 1, :);
        EXk_a = mean(Xk_a);
        EXk_b = mean(Xk_b);
        for i = 1:length(Xk_a)
            cov_Xk(tau, k) = cov_Xk(tau, k) + (Xk_a(i) - EXk_a) * (Xk_b(i) - EXk_b);
        end
        cov_Xk(tau, k) = cov_Xk(tau, k) / length(Xk_a);
    end
end

% vypocet parametru
cov_Xk_vypocet = zeros(TAU, K);
for tau = 1:TAU
    cov_Xk_vypocet(tau, :) = Q*exp((-tau + 1) * b);
end

% Vykresleni
figure
hold on
p1 = plot(0:length(cov_Xk)-1, cov_Xk(1:6, :),'LineWidth', 1.5);
p2 = plot(0:length(cov_Xk_vypocet)-1, cov_Xk_vypocet(1:6, :), 'black:','LineWidth',1.5);
title('Autokovarianèní funkce')
xlabel('Èas [k]')
legend([p1(1) p1(2) p1(3) p1(4) p1(5) p1(6) p2(1)], 'cov[X_k,X_k]', 'cov[X_k,X_{k+1}]', 'cov[X_k,X_{k+2}]', 'cov[X_k,X_{k+3}]', 'cov[X_k,X_{k+4}]', 'cov[X_k,X_{k+5}]', 'Vypoètené')
%% 2
clear
M = 10^4;
N = 100;
TAU = 6; %v zadani indexovano od 0 (6 hodnot)
K = 95; %v zadani indexovano od 0 (95 hodnot)

% generovani realizaci
Xk = zeros(N, M);
for i = 1:M
    for j = 2:N % X0 = 0
        Xk(j, i) = Xk(j - 1, i) + randn; 
    end
end

% vykresleni 8 realizaci
figure
plot(0:N-1, Xk(:, 1:8),'LineWidth', 1.5)
title('8 realizací Wienerova procesu')
xlabel('Èas [k]')
ylabel('X')

% odhad parametru
cov_Xk = zeros(TAU, K);
for tau = 1:TAU
    for k = 1:K
        Xk_a = Xk(k, :);
        Xk_b = Xk(k + tau - 1, :);
        EXk_a = mean(Xk_a);
        EXk_b = mean(Xk_b);
        for i = 1:length(Xk_a)
            cov_Xk(tau, k) = cov_Xk(tau, k) + (Xk_a(i) - EXk_a) * (Xk_b(i) - EXk_b);
        end
        cov_Xk(tau, k) = cov_Xk(tau, k) / length(Xk_a);
    end
end

% vypocet parametru
cov_Xk_vypocet = zeros(TAU, K);
for tau = 1:TAU
    for k = 1:K
        cov_Xk_vypocet(tau, k) = k - 1;
    end
end

% Vykresleni
figure
hold on
p1 = plot(0:length(cov_Xk)-1, cov_Xk(1:6, :),'LineWidth', 1.5);
p2 = plot(0:length(cov_Xk_vypocet)-1, cov_Xk_vypocet(1:6, :), 'black:','LineWidth',1.5);
title('Autokovarianèní funkce')
xlabel('Èas [k]')
legend([p1(1) p1(2) p1(3) p1(4) p1(5) p1(6) p2(1)], 'cov[X_k,X_k]', 'cov[X_k,X_{k+1}]', 'cov[X_k,X_{k+2}]', 'cov[X_k,X_{k+3}]', 'cov[X_k,X_{k+4}]', 'cov[X_k,X_{k+5}]', 'Vypoètené')

%% 3
clear
M = 10^4;
N = 100;
K = 100; %v zadani indexovano od 0 (99 hodnot)

% generovani realizaci
Xk = zeros(N, M);
Zk = zeros(N, M);
for i = 1:M
    Xk(1, i) = 1 + randn * sqrt(5);
    Zk(1, i) = 5 * Xk(1, i) + randn * sqrt(2);
    for j = 2:N
        Xk(j, i) = 0.95 * Xk(j - 1, i) + 0.5 * randn * sqrt(3);
        Zk(j, i) = 5 * Xk(j, i) + randn * sqrt(2);
    end
end

% odhad parametru
E_Zk = zeros(1, K);
E_Xk = zeros(1, K);
var_Zk = zeros(1, K);
var_Xk = zeros(1, K);

for k = 1:K
    a = Xk(k, :);
    b = Zk(k, :);
    Ex = mean(a);
    Ez = mean(b);
    VARx = 0;
    VARz = 0;
    for i = 1:length(a)
        VARx = VARx + (a(i) - Ex)^2;
        VARz = VARz + (b(i) - Ez)^2;
    end
    VARx = VARx/length(a);
    VARz = VARz/length(b);
    
    E_Xk(k) = Ex;
    var_Xk(k) = VARx;
    E_Zk(k) = Ez;
    var_Zk(k) = VARz;
end

% vypocet parametru
E_Xk_vypocet = zeros(1, K);
var_Xk_vypocet = zeros(1, K);
E_Zk_vypocet = zeros(1, K);
var_Zk_vypocet = zeros(1, K);

for k = 1:K
    E_Xk_vypocet(k) = 0.95^(k-1);
    E_Zk_vypocet(k) = 5 * (0.95^(k-1));
    
    suma = 0;
    for n = 0:k-1
        suma = suma + 0.95^(2*n);
    end
    var_Xk_vypocet(k) = (0.95^(2*(k-1))) * 4 + (0.5^2) * 3 * suma;
    var_Zk_vypocet(k) = var_Xk_vypocet(k) * 5^2 + 2;
end

% vykresleni
figure
hold on
p1 = plot(0:length(E_Xk)-1, E_Xk,'LineWidth',1.5);
p2 = plot(0:length(E_Xk_vypocet)-1, E_Xk_vypocet, 'black:','LineWidth',1.5);
title('Støedni hodnoty X_k')
xlabel('Èas [k]')
ylabel('E[X_k]')
legend([p1(1) p2(1)], 'Odhadované', 'Vypoètené');

figure
hold on
p1 = plot(0:length(E_Zk)-1, E_Zk,'LineWidth',1.5);
p2 = plot(0:length(E_Zk_vypocet)-1, E_Zk_vypocet, 'black:','LineWidth',1.5);
title('Støedni hodnoty Z_k')
xlabel('Èas [k]')
ylabel('E[Z_k]')
legend([p1(1) p2(1)], 'Odhadované', 'Vypoètené')

figure
hold on
p1 = plot(0:length(var_Xk)-1, var_Xk,'LineWidth',1.5);
p2 = plot(0:length(var_Xk_vypocet)-1, var_Xk_vypocet, 'black:','LineWidth',1.5);
title('Variance X_k')
xlabel('Èas [k]')
ylabel('var[X_k]')
legend([p1(1) p2(1)], 'Odhadované', 'Vypoètené')

figure
hold on
p1 = plot(0:length(var_Zk)-1, var_Zk,'LineWidth',1.5);
p2 = plot(0:length(var_Zk_vypocet)-1, var_Zk_vypocet, 'black:','LineWidth',1.5);
title('Variance Z_k')
xlabel('Èas [k]')
ylabel('var[Z_k]')
legend([p1(1) p2(1)], 'Odhadované', 'Vypoètené')


