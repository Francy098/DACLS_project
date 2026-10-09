function [E_f, Z_f, Hw_optF, E_b, Z_b, Hw_optB, E_c, Z_c, Hw_optC] = Optimal_Filtering(tau, L, Y1, V_noise)
% OPTIMAL_FILTERING Calcola i filtri ottimi (Wiener / Least Square in forma chiusa)
% per le configurazioni Linear/Forward, Backward e Combined Approach.

[K, M] = size(Y1); 
alpha = 0.95;

start_m_f = 0 + tau + L;

%% 1. Approccio Forward Least Square
Hw_optF = zeros(L, K, M);
Ryy_f   = zeros(L, L, K, M); 
ryv_f   = zeros(L, K, M); 
E_f     = Y1; 
Z_f     = zeros(K, M);

for k = 1:K
    for m = start_m_f:M
        % Vettore finestra regressiva
        Yvec = Y1(k, (m-tau):-1:(m-tau-L+1)).';
        
        % Matrice di autocorrelazione Ryy
        Ryy_f(:, :, k, m) = alpha * Ryy_f(:, :, k, m-1) + (1-alpha) * (Yvec * Yvec');
        reg = 1e-6 * trace(Ryy_f(:, :, k, m)) / L + 1e-12;  % Per pesare i valori su tutti i campioni
        Ryy_temp = Ryy_f(:, :, k, m) + reg * eye(L);
        
        % Cross-correlazione
        ryv_f(:, k, m) = alpha * ryv_f(:, k, m-1) + (1-alpha) * Yvec * conj(V_noise(k, m)); 
        
        % Filtro ottimo di Wiener
        Hw_optF(:, k, m) = Ryy_temp \ ryv_f(:, k, m);
        
        Z_f(k, m) = Hw_optF(:, k, m)' * Yvec;
        E_f(k, m) = Y1(k, m) - Z_f(k, m);
    end
end

%% 2. Approccio Backward Least Square (BLS)
Hw_optB = zeros(L, K, M); 
Ryy_b   = zeros(L, L, K, M+1); 
ryv_b   = zeros(L, K, M+1); 
E_b     = Y1; 
Z_b     = zeros(K, M);

start_m_b = M - tau - L + 1;

for k = 1:K
    for m = start_m_b:-1:1
        % Vettore finestra progressiva di Y1 (L elementi da m+tau a m+tau+L-1)
        Yvec = Y1(k, (m+tau):1:(m+tau+L-1)).';
        
        % Matrice di autocorrelazione Ryy (propagazione all'indietro m+1)
        Ryy_b(:, :, k, m) = alpha * Ryy_b(:, :, k, m+1) + (1-alpha) * (Yvec * Yvec');
        reg = 1e-6 * trace(Ryy_b(:, :, k, m)) / L + 1e-12;  
        Ryy_temp = Ryy_b(:, :, k, m) + reg * eye(L);
        
        % Cross-correlazione
        ryv_b(:, k, m) = alpha * ryv_b(:, k, m+1) + (1-alpha) * Yvec * conj(V_noise(k, m)); 
        
        % Filtro ottimo di Wiener
        Hw_optB(:, k, m) = Ryy_temp \ ryv_b(:, k, m);
        
        Z_b(k, m) = Hw_optB(:, k, m)' * Yvec;
        E_b(k, m) = Y1(k, m) - Z_b(k, m);
    end
end

%% 3. Approccio Combinato -> Selezione a Minima Energia
Hw_optC = zeros(L, K, M); 
E_c     = zeros(K, M);
Z_c     = zeros(K, M);

en_Eb = abs(E_b).^2;
en_Ef = abs(E_f).^2;
en_Y1 = abs(Y1).^2;

for k = 1:K
    for m = 1:M
        if (en_Eb(k, m) <= en_Ef(k, m)) && (en_Eb(k, m) <= en_Y1(k, m))
            E_c(k, m) = E_b(k, m); 
            Z_c(k, m) = Z_b(k, m);
            Hw_optC(:, k, m) = Hw_optB(:, k, m);

        elseif (en_Ef(k, m) <= en_Eb(k, m)) && (en_Ef(k, m) <= en_Y1(k, m))
            E_c(k, m) = E_f(k, m); 
            Z_c(k, m) = Z_f(k, m);
            Hw_optC(:, k, m) = Hw_optF(:, k, m);

        else
            E_c(k, m) = Y1(k, m); 
            Z_c(k, m) = 0;
            Hw_optC(:, k, m) = zeros(L, 1);
        end
    end
end

end