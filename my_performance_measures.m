function [zita_nr_full, vsd_full, NRR] = my_performance_measures(tau, L, Hw, V_noise, X_clean, L_bank)

[K, M] = size(X_clean(:,1:449));
alpha = 0.95;

% Preallocazione varianze e cross-correlazioni principali
phi_X = zeros(K, M);
phi_V = zeros(K, M);

% Matrici 3D per le covarianze dell'interferenza [L x L x K]
Phi_Xint = zeros(L, L, K);
Phi_Vint = zeros(L, L, K);

% Preallocazione segnali e distorsioni
Xfd = zeros(K, M);
Xri = zeros(K, M);
Vrn = zeros(K, M);
phi_Xfd = zeros(K, M);
phi_Xri = zeros(K, M);
phi_Vrn = zeros(K, M);

% Cross-correlazioni ricorsive [K x L]
cross_x = zeros(K, L);
cross_v = zeros(K, L);

% Vettori per accumulare i valori scalari medi di gamma
gamma_x = zeros(L,K,M);
gamma_v = zeros(L,K,M);

% Preallocazione metriche a NaN per ignorare il transitorio nella media
vsd_km = zeros(K, M);
zita_nr = NaN(K, M);
den_zita = NaN(K, M); 

%% 1. CALCOLO FRAME-BY-FRAME ED EQUAZIONI (28)-(46)
for m = (tau + L):M
    for k = 1:K
        l_curr = L_bank(k, m);
        
        % Controllo di sicurezza per filtri disattivati o di lunghezza nulla
        if l_curr <= 0
            continue;
        end
        
        % --- A. Stima varianze temporali phi_X e phi_V ---
        phi_X(k, m) = alpha * phi_X(k, m-1) + (1-alpha) * abs(X_clean(k, m))^2;
        phi_V(k, m) = alpha * phi_V(k, m-1) + (1-alpha) * abs(V_noise(k, m))^2;
        
        % --- B. Costruzione dei vettori gamma_x e gamma_v per ritardi [tau ... tau + l_curr - 1] (Eq. 13) ---
        gamma_x_vec = zeros(l_curr, 1);
        gamma_v_vec = zeros(l_curr, 1);
        X_int_vec   = zeros(l_curr, 1);
        V_int_vec   = zeros(l_curr, 1);
        
        for r = 0:(l_curr - 1)
            delay_idx = m - (tau + r) + 1;
            
            if delay_idx >= 1
                % Cross-correlazione ricorsiva coerente sull'indice r+1
                cross_x(k, r+1) = alpha * cross_x(k, r+1) + (1-alpha) * X_clean(k, m) * conj(X_clean(k, delay_idx));
                cross_v(k, r+1) = alpha * cross_v(k, r+1) + (1-alpha) * V_noise(k, m) * conj(V_noise(k, delay_idx));
                
                % Normalizzazione rispetto alla varianza del frame corrente m (stabilità numerica)
                gamma_x_vec(r+1) = cross_x(k, r+1) / (phi_X(k, m) + 1e-6);
                gamma_v_vec(r+1) = cross_v(k, r+1) / (phi_V(k, m) + 1e-6);               
                
                % Limitazione del modulo a 1 per prevenire la crescita esponenziale
                if abs(gamma_x_vec(r+1)) > 1
                    gamma_x_vec(r+1) = gamma_x_vec(r+1) / abs(gamma_x_vec(r+1));
                end
                if abs(gamma_v_vec(r+1)) > 1
                    gamma_v_vec(r+1) = gamma_v_vec(r+1) / abs(gamma_v_vec(r+1));
                end
                % Componenti interferenti x_int e v_int
                X_int_vec(r+1) = X_clean(k, delay_idx) - conj(gamma_x_vec(r+1)) * X_clean(k, m);
                V_int_vec(r+1) = V_noise(k, delay_idx) - conj(gamma_v_vec(r+1)) * V_noise(k, m);
            end
        end
        
        % Accumulo scalare per output
        gamma_x(1:l_curr,k,m) = abs(gamma_x_vec);
        gamma_v(1:l_curr,k,m) = abs(gamma_v_vec);
        
        % --- C. Aggiornamento ricorsivo matrici di covarianza Phi_X' e Phi_V' (Eq. 29-30) ---
        Phi_Xint(1:l_curr, 1:l_curr, k) = alpha * Phi_Xint(1:l_curr, 1:l_curr, k) + (1-alpha) * (X_int_vec * X_int_vec');
        Phi_Vint(1:l_curr, 1:l_curr, k) = alpha * Phi_Vint(1:l_curr, 1:l_curr, k) + (1-alpha) * (V_int_vec * V_int_vec');
        
        % --- D. Estrazione filtro attivo h di dimensione l_curr ---
        h = Hw(1:l_curr, k, m);
        
        % --- E. Componenti dei segnali filtrati e residui ---
        Xfd(k, m) = X_clean(k, m) * (1 - h' * gamma_x_vec);
        Xri(k, m) = -h' * X_int_vec;
        Vrn(k, m) = V_noise(k, m) * (1 - h' * gamma_v_vec) - (h' * V_int_vec);
        
        % --- F. Potenze teoriche e Distorsione ---
        tmp_x = 1 - h' * gamma_x_vec;
        phi_Xfd(k, m) = phi_X(k, m) * abs(tmp_x)^2;
        
        phi_Xri(k, m) = real(h' * Phi_Xint(1:l_curr, 1:l_curr, k) * h);
        
        tmp_v = 1 - h' * gamma_v_vec;
        phi_Vrn(k, m) = phi_V(k, m) * abs(tmp_v)^2 + real(h' * Phi_Vint(1:l_curr, 1:l_curr, k) * h);
        
        %Calcolo vsd (Formula 65)
        vsd_km(k,m) = abs(h' * conj(gamma_x_vec))^2;

    end
end

%% 2. RIDUZIONE A SCALARE ED ENERGIA GLOBALE FULL-BAND

% Metriche Full-Band
start_frame = tau + L;

vsd_full = mean(mean(vsd_km(:, start_frame:end), 1, 'omitnan'), 'omitnan');

% Calcolo zita_nr Formula (40)
den_zita2 = sum(phi_Xri, 1,'omitnan') + sum(phi_Vrn, 1,'omitnan');
num_zita2 = sum(phi_V, 1,'omitnan');
zita_nr2_temp = num_zita2 ./ den_zita2;
zita_nr_full = mean(zita_nr2_temp, 'omitnan');

%% 3. CALCOLO NRR (Noise Rejection Ratio) 
NRR = zeros(1, M);
E_v_in  = sum(abs(V_noise).^2, 1); % Energia rumore in ingresso
E_v_out = sum(abs(Vrn).^2, 1);     % Energia residua in uscita

% NRR = 10 * log10 ( Potenza IN / Potenza OUT )
NRR = 10 * log10((E_v_in + 1e-12) ./ (E_v_out + 1e-12)); % +50dB, aggiustamento offset

end