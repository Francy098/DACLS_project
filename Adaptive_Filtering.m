function [E_f,Z_f,Hw_f, E_b,Z_b,Hw_b, E_c,Z_c,Hw_c, Lf_bank,Lb_bank,Lc_bank] = Adaptive_Filtering(tau,L,Y1_temp,V_noise)

Y1 = Y1_temp(:,1:449); 
[K,M] = size(Y1);

% Inizializzazione parametri
u = 0.5 * ones(K,1); % Step size per bin di frequenza
delta = 1e-4;

%% 1. Noise Detector (Vettorializzato)
V_max = max(mag2db(abs(V_noise) + 1e-6), [], 'all');
Itresh = -25; % Soglia in dB
int = (mag2db(abs(V_noise) + 1e-6) - V_max) > Itresh; % Matrice logica K  x  M

%% 2. Approccio Forward FMLNLMS
Hw_f = zeros(L, K, M+1);
E_f = Y1; % Default trasparente
Z_f = zeros(K, M);
Lf_bank = ones(K, M);

for k = 1:K
    % Garantisce che m - tau - L + 1 >= 1
    start_m = max(tau + L, tau + 1); 
    for m = start_m:(M-1)
        % Mapped-L Forward
        l_f = 1;
        for i = 0:L-1
            idx_check = m - tau - i;
            if idx_check >= 1 && ...
                    int(k, idx_check)
                l_f = l_f + 1;
            else
                break;
            end
        end
        l_f = min(l_f, L); % l_f non può mai superare L
        Lf_bank(k,m) = l_f;
        
        % Vettore regressore sicuro
        idx_end = m - tau + 1;
        idx_start = idx_end - l_f + 1;
        Y2 = Y1(k, idx_end:-1:idx_start).';
        
        Z_f(k,m) = Hw_f(1:l_f, k, m)' * Y2;
        E_f(k,m) = Y1(k,m) - int(k,m) * Z_f(k,m);
        
        den_NLMS = Y2'*Y2 + delta;
        if int(k,m) % ==1 allora è presente il rumore
            Hw_f(1:l_f, k, m+1) = Hw_f(1:l_f, k, m) + (u(k) * Y2 * conj(E_f(k,m))) / den_NLMS;
        else
            Hw_f(1:l_f, k, m+1) = Hw_f(1:l_f, k, m);
        end
    end
end

%% 3. Approccio Backward BMLNLMS
Hw_b = zeros(L, K, M+1);
E_b = Y1; 
Z_b = zeros(K, M);
Lb_bank = ones(K, M);

for k = 1:K
    % Garantisce che m + tau + L - 1 <= M
    start_m = min(M - tau, M - tau - L + 1);
    for m = start_m:-1:1
        % Mapped-L Backward
        l_b = 1;
        for i = 0:L-1
            idx_check = m + tau + i;
            if idx_check <= M && ...
                    int(k, idx_check)
                l_b = l_b + 1;
            else
                break;
            end
        end
        l_b = min(l_b, L);
        Lb_bank(k,m) = l_b;
        
        % Vettore regressore sicuro
        idx_start = m + tau;
        idx_end = min(idx_start + l_b - 1, M);
        Y2 = Y1(k, idx_start:idx_end).';
        actual_l = length(Y2); % Controllo dimensioni effettive
        % backward recursion: h(m) viene ottenuto da h(m+1)
        Z_b(k,m) = Hw_b(1:actual_l, k, m+1)' * Y2; 
        E_b(k,m) = Y1(k,m) - int(k,m) * Z_b(k,m);
        
        den_NLMS = Y2'*Y2 + delta;
        if int(k,m)
            Hw_b(1:actual_l, k, m) = Hw_b(1:actual_l, k, m+1) + (u(k) * Y2 * conj(E_b(k,m))) / den_NLMS;
        else
            Hw_b(1:actual_l, k, m) = Hw_b(1:actual_l, k, m+1);
        end
    end
end

%% 4. Approccio Combinato -> selezione a Minima Energia

%Inizializzazione
Hw_c = zeros(L,K,M+1);
E_c = zeros(K,M);
Z_c = zeros(K,M);
Lc_bank = ones(K,M);

% Energie spettrali
en_Eb = zeros(K,M);
en_Ef = zeros(K,M);
en_Y1 = zeros(K,M);

for k = 1:K %Iterazione sulle RIGHE

    for m = 1:M
        en_Eb(k,m) = abs(E_b(k,m)).^2;
        en_Ef(k,m) = abs(E_f(k,m)).^2;
        en_Y1(k,m) = abs(Y1(k,m)).^2;
            
      
        if((en_Eb(k,m) <= en_Ef(k,m)) && ((en_Eb(k,m) <= en_Y1(k,m))))
            E_c(k,m) = E_b(k,m); %Backward
            Z_c(k,m) = Z_b(k,m);
            
            Hw_c(:,k,m) = Hw_b(:,k,m);
            Lc_bank(k,m) = Lb_bank(k,m);

        elseif((en_Eb(k,m) > en_Ef(k,m)) && ((en_Ef(k,m) <= en_Y1(k,m))))
            E_c(k,m) = E_f(k,m); %Forward
            Z_c(k,m) = Z_f(k,m);
            
            Hw_c(:,k,m) = Hw_f(:,k,m);
            Lc_bank(k,m) = Lf_bank(k,m);

        else
            E_c(k,m) = Y1(k,m); % No noise 
            Z_c(k,m) = 0;

            Hw_c(:,k,m) = zeros(L,1);
            Lc_bank(k,m) = 0;
        end

        Z_c(k,m) = Y1(k,m) - E_c(k,m); % ci serve per fare un'analisi delle prestazioni
    end
end

end