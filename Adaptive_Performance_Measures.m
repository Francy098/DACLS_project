function [zita_nr_full, vsd_full, NRR] = Adaptive_Performance_Measures(tau, L, Hw, V_noise, X_clean, L_bank)

[K, M] = size(X_clean(:,1:449));
alpha = 0.95;

% Preallocazione varianze e cross-correlazioni principali
phi_X = zeros(K, M);
phi_V = zeros(K, M);

% Preallocazione segnali e distorsioni
Xfd = zeros(K, M);
Xri = zeros(K, M);
Vrn = zeros(K, M);
phi_Xfd = zeros(K, M);
phi_Xri = zeros(K, M);
phi_Vrn = zeros(K, M);

% Preallocazione metriche a NaN per ignorare il transitorio nella media
vsd_km = zeros(K, M); 

%% 1. CALCOLO FRAME-BY-FRAME ED EQUAZIONI (28)-(46)
start_frame = tau + L;

for k = 1:K

    PHI_X1 = zeros(tau+L, tau+L);
    PHI_V1 = zeros(tau+L, tau+L);

    for m = (tau + L):M

        l_curr = L_bank(k, m);
        
        % Controllo di sicurezza per filtri disattivati o di lunghezza nulla
        if l_curr <= 0
            continue;
        end
        
        % --- A. Stima varianze temporali phi_X e phi_V ---
        phi_X(k, m) = sum(abs(X_clean(k, m : -1 : m-l_curr+1)).^2) / l_curr;
        phi_V(k, m) = sum(abs(V_noise(k, m : -1 : m-l_curr+1)).^2) / l_curr;
        
        % Windowing X_clean
        X1 = X_clean(k, m : -1 : m-tau-L+1); % X1 is a Window of X_clean
        X1 = X1(:); % Vettore colonna
        %PHI_X1 = xcorr(X1, X1','biased');
        PHI_X1 = alpha * PHI_X1 + (1-alpha) * (X1 * X1') / L;
        % Signal correlation vector (1xL)
        gamma_tmp = PHI_X1(1, tau+1 : (tau+l_curr)) / (PHI_X1(1,1) + 1e-12); 
        gamma_x = gamma_tmp(:); % ci assicura che sia un vettore colonna
        PHI_X_dly =  PHI_X1(tau+1 : tau+l_curr, tau+1 : tau+l_curr); %delayed correlation matrix

        % Windowing V_noise
        V1 = V_noise(k, m:-1:m-tau-L+1); % X1 is a Window of X_clean
        V1 = V1(:); % Vettore colonna
        %PHI_V1 = xcorr(V1, V1','biased');
        PHI_V1 = alpha * PHI_V1 + (1-alpha) * (V1 * V1') / L;
        % Noise correlation vector (1xL)
        gamma_tmp = PHI_V1(1, tau+1 : (tau+l_curr)) / (PHI_V1(1,1)+ 1e-12); 
        gamma_v = gamma_tmp(:); % ci assicura che sia un vettore colonna
        PHI_V_dly =  PHI_V1(tau+1 : tau+l_curr, tau+1 : tau+l_curr); %delayed correlation matrix
        
        % Signal interference vector (1xL)
        X_clean_vec = X_clean(k, m-tau :-1: m-tau-l_curr+1);
        X_clean_vec = X_clean_vec(:); 
        X_int_vec = X_clean_vec - conj(gamma_x) * X_clean(k,m); % garantiamo sia vettore colonna
        % Noise interference vector (1xL)
        V_noise_vec = V_noise(k, m-tau :-1: m-tau-l_curr+1);
        V_noise_vec = V_noise_vec(:);
        V_int_vec = V_noise_vec - conj(gamma_v) * V_noise(k,m); % garantiamo sia vettore colonna
        
        % --- C. Matrici di covarianza Phi_X' e Phi_V' (Eq. 33-34) ---
        Phi_Xint = PHI_X_dly - phi_X(k,m) * conj(gamma_x) * conj(gamma_x)';
        Phi_Vint = PHI_V_dly - phi_V(k,m) * conj(gamma_v) * conj(gamma_v)';
        
        % --- D. Estrazione filtro attivo h di dimensione l_curr ---
        h = Hw(1:l_curr, k, m);
        
        % --- E. Componenti dei segnali filtrati e residui ---
        Xfd(k, m) = X_clean(k, m) * (1 - h' * conj(gamma_x));
        Xri(k, m) = -h' * X_int_vec;
        Vrn(k, m) = V_noise(k, m) * (1 - h' * conj(gamma_v)) - (h' * V_int_vec);
        
        % --- F. Potenze teoriche e Distorsione ---
        tmp_x = 1 - h' * conj(gamma_x);
        phi_Xfd(k, m) = phi_X(k, m) * abs(tmp_x)^2;
        
        phi_Xri(k, m) = h' * Phi_Xint * h;
        
        tmp_v = 1 - h' * conj(gamma_v);
        phi_Vrn(k, m) = phi_V(k, m) * abs(tmp_v)^2 + h' * Phi_Vint * h;
        
        %Calcolo metriche
        vsd_km(k,m) = abs(h' * conj(gamma_x))^2;
        %vsd_num(k, m) = alpha * vsd_num(k, m-1) + (1-alpha) * abs(X_clean(k, m) - Xfd(k, m))^2;
        
        %zita_pt1 = real(h' * Phi_Xint(1:l_curr,1:l_curr,k)*h); 
        %zita_pt2 = phi_V(k,m) * abs(1 - h' * gamma_v_vec)^2;
        %zita_pt3 = real(h' * Phi_Vint(1:l_curr,1:l_curr,k) * h);
        %den_zita = zita_pt1 + zita_pt2 + zita_pt3 + 1e-12;
        %zita_nr(k,m) = phi_V(k,m) * 1/den_zita;

    end
end

%% 2. RIDUZIONE A SCALARE ED ENERGIA GLOBALE FULL-BAND

% =========================================================================
% VECCHIA IMPLEMENTAZIONE (COMMENTATA)
% =========================================================================
% vsd_temp = sum(vsd_num, 1, 'omitnan') ./ (sum(phi_X, 1, 'omitnan') + 1e-12);
% for i = 1 : 449
%     if vsd_temp(i) > = 1
%         vsd_temp(i) = NaN;
%     end
% end
% vsd_full = mean(vsd_temp, 'omitnan');
%
% den_zita2 = sum(phi_Xri, 1,'omitnan') + sum(phi_Vrn, 1,'omitnan');
% num_zita2 = sum(phi_V, 1,'omitnan');
% zita_nr2_temp = num_zita2 ./ den_zita2;
% zita_nr_full = mean(zita_nr2_temp, 'omitnan');
%
% E_v_in  = sum(phi_V, 1, 'omitnan');
% E_v_out = sum(phi_Vrn, 1, 'omitnan');
% NRR = 10 * log10((E_v_in + 1e-12) ./ (E_v_out + 1e-12));
% =========================================================================

% Metriche Full-Band
vsd_full = mean(mean(vsd_km(:, start_frame:end), 1, 'omitnan'), 'omitnan');

% Calcolo zita_nr Formula (40)
den_zita2 = sum(phi_Xri, 1,'omitnan') + sum(phi_Vrn, 1,'omitnan');
num_zita2 = sum(phi_V, 1,'omitnan');
zita_nr2_temp = num_zita2 ./ den_zita2;
zita_nr_full = mean(zita_nr2_temp, 'omitnan');

%% 3. CALCOLO NRR (Noise Rejection Ratio) 
NRR = zeros(1, M);
% ----------------------------------------------------------------------------------
% Prima era (start_frame:M) perché cercavamo di evitare il transittorio
%
E_v_in  = sum(abs(V_noise).^2, 1); % Potenza rumore in ingresso
E_v_out = sum(abs(Vrn).^2, 1);     % Potenza residua in uscita

% NRR = 10 * log10 ( Potenza IN / Potenza OUT )
NRR = 10 * log10((E_v_in + 1e-12) ./ (E_v_out + 1e-12)) + 50; % +50dB, aggiustamento offset
end