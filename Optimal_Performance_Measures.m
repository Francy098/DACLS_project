function [zita_nr_full, vsd_full, NRR] = Optimal_Performance_Measures(tau, L, Hw, V_noise, X_clean)

[K, M] = size(X_clean);
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
vsd_km = NaN(K, M);
phi_Vrn = NaN(K, M);
Vrn = zeros(K, M);

%% 1. CALCOLO FRAME-BY-FRAME ED EQUAZIONI (28)-(46)
start_frame = tau + L;

for k = 1:K

    PHI_X1 = zeros(tau+L, tau+L);
    PHI_V1 = zeros(tau+L, tau+L);

    for m = (start_frame):M

        
        % --- A. Stima varianze temporali phi_X e phi_V ---
        %phi_X(k, m) = alpha * phi_X(k, m-1) + (1-alpha) * abs(X_clean(k, m))^2;
        %phi_V(k, m) = alpha * phi_V(k, m-1) + (1-alpha) * abs(V_noise(k, m))^2;

        % Implementazione Estimed autocorrelation del prof
        phi_X(k, m) = sum(abs(X_clean(k, m : -1 : m-L+1)).^2) / L;
        phi_V(k, m) = sum(abs(V_noise(k, m : -1 : m-L+1)).^2) / L;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inizializzazione degli accumulatori
        %temp_phiX = 0;
        %temp_phiV = 0;
        % Accumulo del modulo al quadrato lungo la finestra L
        %for i = 0 : L-1
        %    temp_phiX = temp_phiX + abs(X_clean(k, m - i))^2;
        %    temp_phiV = temp_phiV + abs(V_noise(k, m - i))^2;
        %end
        % Calcolo della media temporale su L campioni (Slide 8)
        %phi_X(k, m) = temp_phiX / L;
        %phi_V(k, m) = temp_phiV / L;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Windowing X_clean
        X1 = X_clean(k, m : -1 : m-tau-L+1); % X1 is a Window of X_clean
        %X1 = X1(:); % Vettore colonna
        
        % Calculating PHI_X1 correlation matrix
        %PHI_X1 = alpha * PHI_X1 + (1-alpha) * (X1 * X1') / L;
        PHI_vec = xcorr(X1, X1','biased');
        num_zeros = L+tau- length(PHI_vec);
        PHI_X1_temp = [xcorr(X1, X1','biased'), zeros(1,num_zeros)];
        idx_start = 1;
        idx_end = L + tau;

        for i = 1 : L+tau
            PHI_X1(i, 1 : (L+tau)) = PHI_X1_temp(idx_start : idx_end);
            idx_start = idx_start + (L + tau);
            idx_end = idx_end + (L + tau);
        end
        
        % Signal correlation vector (1xL)
        gamma_tmp = PHI_X1(1, tau+1 : (tau+L)) / (PHI_X1(1,1) + 1e-12); % -1 perché deve essere lungo quanto L!!
        gamma_x = gamma_tmp(:); % ci assicura che sia un vettore colonna
        PHI_X_dly =  PHI_X1(tau+1 : tau+L, tau+1 : tau+L); %delayed correlation matrix

        % Windowing V_noise
        V1 = V_noise(k, m:-1:m-tau-L+1); % X1 is a Window of X_clean
        %V1 = V1(:); % Vettore colonna
        %PHI_V1 = alpha * PHI_V1 + (1-alpha) * (V1 * V1') / L;
        PHI_vec = xcorr(V1, V1','biased');
        num_zeros = L+tau- length(PHI_vec);
        PHI_V1_temp = [xcorr(V1, V1','biased'), zeros(1,num_zeros)];

        idx_start = 1;
        idx_end = L + tau;

        for i = 1 : L+tau
            PHI_V1(i, 1 : (L+tau)) = PHI_V1_temp(idx_start : idx_end);
            idx_start = idx_start + (L + tau);
            idx_end = idx_end + (L + tau);
        end
        % Noise correlation vector (1xL)
        gamma_tmp = PHI_V1(1, tau+1 : (tau+L)) / (PHI_V1(1,1)+ 1e-12); % vettore (1x(tau+L-1)) % -1 perché deve essere lungo quanto L!!
        gamma_v = gamma_tmp(:); % ci assicura che sia un vettore colonna
        PHI_V_dly =  PHI_V1(tau+1 : tau+L, tau+1 : tau+L); %delayed correlation matrix
        
        % Signal interference vector (1xL)
        X_clean_vec = X_clean(k, m-tau :-1: m-tau-L+1);
        X_clean_vec = X_clean_vec(:); 
        X_int_vec = X_clean_vec - conj(gamma_x) * X_clean(k,m); % garantiamo sia vettore colonna
        % Noise interference vector (1xL)
        V_noise_vec = V_noise(k, m-tau :-1: m-tau-L+1);
        V_noise_vec = V_noise_vec(:);
        V_int_vec = V_noise_vec - conj(gamma_v) * V_noise(k,m); % garantiamo sia vettore colonna
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% LEGENDA
        %%% Dove vediamo scritto sul paper 
        % - (k,m-tau) ---> vettore di 1xL elementi da(m-tau : -1 : m-tau-L+1) 
        % - (k,m,tau) ---> matrice LxL che VARIA OGNI (k,m)     
        % - conj(--)' ---> trasposta pur

        % --- C. Estrazione filtro attivo h di dimensione L ---
        h = Hw(:, k, m);
        
        % --- D. Componenti dei segnali filtrati e residui ---
        Xfd(k, m) = X_clean(k, m) * (1 - h' * conj(gamma_x)); % ----------- inserito il CONJ!!!!!!!
        Xri(k, m) = -h' * X_int_vec;
        Vrn(k, m) = V_noise(k, m) * (1 - h' * conj(gamma_v)) - (h' * V_int_vec); % ----------- inserito il CONJ!!!!!!!
        
        % --- E. Aggiornamento ricorsivo matrici di covarianza Phi_X' e Phi_V' (Eq. 29-30) ---
        %Phi_Xint(:, :, k) = alpha * Phi_Xint(:, :, k) + (1-alpha) * (X_int_vec * X_int_vec');
        %PHI_Vint(:, :, k) = alpha * PHI_Vint(:, :, k) + (1-alpha) * (V_int_vec * V_int_vec');
        % ---- Nuovo metodo ----
        %Phi_Xint(:, :, k) = xcorr(X_int_vec, X_int_vec', 'biased'); % cross-correlazione mediata sulla lunghezza dei vettori
        %PHI_Vint(:, :, k) = xcorr(V_int_vec, V_int_vec', 'biased'); % cross-correlazione mediata sulla lunghezza dei vettori

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Provare il calcolo di Phi_Xint e Phi_Vint come indicato nella
        % (33) e (34) ??
        Phi_Xint = PHI_X_dly - phi_X(k,m) * conj(gamma_x) * conj(gamma_x)';
        Phi_Vint = PHI_V_dly - phi_V(k,m) * conj(gamma_v) * conj(gamma_v)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % --- F. Potenze teoriche e Distorsione ---
        tmp_x = 1 - h' * conj(gamma_x); % ----------- inserito il CONJ!!!!!!!
        phi_Xfd(k, m) = phi_X(k, m) * abs(tmp_x)^2;
        
        phi_Xri(k, m) = h' * Phi_Xint * h; % ----------- tolto il REAL!!!!!!!
        
        tmp_v = 1 - h' * conj(gamma_v);         % ----------- inserito il CONJ!!!!!!!
        phi_Vrn(k, m) = phi_V(k, m) * abs(tmp_v)^2 + h' * Phi_Vint * h;
        
        %Calcolo vsd (Formula 45)
        vsd_km(k, m) = abs(h' * conj(gamma_x))^2;
        
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
E_v_in  = sum(abs(V_noise).^2, 1); % Potenza rumore in ingresso
E_v_out = sum(abs(Vrn).^2, 1);     % Potenza residua in uscita

% NRR = 10 * log10 ( Potenza IN / Potenza OUT ) +50dB, aggiustamento offset
NRR = 10 * log10((E_v_in + 1e-12) ./ (E_v_out + 1e-12)) + 50; 

end