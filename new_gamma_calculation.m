function [gamma_x_mean, gamma_v_mean] = new_gamma_calculation(K, M, L, Fs)

%% Estrazione file audio dal DATABASE
%
% Qui calcoliamo il gamma_x mediato su tutti i 20 speech nel database,
% mentre il gamma_v lo mediamo su tutti i 26 rumori nel database.
%
[all_speeches, all_noises] = my_all_database(Fs);
%
num_speeches = size(all_speeches, 2); 
num_tau = 26; % tau da 0 a 25

alpha = 0.95;
gamma_x_all = zeros(num_tau, num_speeches);

% Parametri STFT
win_length = 512;
win = hamming(win_length); 
overlap = round(win_length * 0.75); 
nfft = 512;


%% 1. CALCOLO GAMMA_X SUL DATABASE

for speech_idx = 1:num_speeches 

    x_clean = all_speeches(:, speech_idx);
    X_clean = stft(x_clean, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap); 
    
    
    for tau_idx = 1:num_tau
        tau = tau_idx - 1; % tau va da 0 a 25
        
        % Reset variabili per ogni nuova combinazione di speech e tau
        phi_X = zeros(K, M);
        gamma_x_mat = NaN(L, K, M); % Preallocazione a NaN per ignorare i frame non calcolati
        
        start_frame = tau + L;
        
        for k = 1:K
        
            %PHI_X1 = zeros(tau+L, tau+L);
        
            for m = (start_frame):M
        
                % Implementazione Estimed autocorrelation del prof
                phi_X(k, m) = sum(abs(X_clean(k, m : -1 : m-L+1)).^2) / L;
                
                % Windowing X_clean
                X1 = X_clean(k, m : -1 : m-tau-L+1); % X1 is a Window of X_clean
                %X1 = X1(:); % Vettore colonna
                %PHI_X1 = alpha * PHI_X1 + (1-alpha) * (X1 * X1') / L;
                PHI_X1 = xcorr(X1, X1','biased');
                % Signal correlation vector (1xL)
                gamma_tmp = PHI_X1(1, tau+1 : (tau+L)) / (PHI_X1(1,1) + 1e-12); % -1 perché deve essere lungo quanto L!!
                gamma_x = gamma_tmp(:); % ci assicura che sia un vettore colonna
                
                % Salvataggio gamma_x corrente
                gamma_x_mat(:, k, m) = abs(gamma_x);
            end
        end
           % Calcolo media scalare ignorando i NaN del transitorio iniziale
           gamma_x_all(tau_idx, speech_idx) = mean(gamma_x_mat, 'all', 'omitnan');
    end
end

%% 2. CALCOLO GAMMA_V SUL DATABASE RUMORI


num_noises = size(all_noises, 2); % numero di rumori 26

gamma_v_all = zeros(num_tau, num_noises);

for noise_idx = 1:num_noises 
    
    v_noise = all_noises(:, noise_idx);
    V_noise = stft(v_noise, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap); 
    [K_noise, M_noise] = size(V_noise);           % Estrazione dimensioni reali
    

    
    for tau_idx = 1:num_tau
        tau = tau_idx - 1;
        
        % Allocazione corretta con le dimensioni reali di V_noise
        phi_V = zeros(K_noise, M_noise);
        gamma_v_mat = NaN(L, K_noise, M_noise);
        
        start_frame = tau + L;

        for k = 1:K_noise
        
            %PHI_V1 = zeros(tau+L, tau+L);

            % Inizializzazione della matrice al primo frame per evitare divisioni per zero
            %V1_init = V_noise(k, start_frame : -1 : start_frame - tau - L + 1);
            %V1_init = V1_init(:);
            %PHI_V1 = (V1_init * V1_init') / L;
        
            for m = (start_frame) : M_noise
        
                % Implementazione Estimed autocorrelation del prof
                phi_V(k, m) = sum(abs(V_noise(k, m : -1 : m-L+1)).^2) / L;

                % Windowing V_noise
                V1 = V_noise(k, m:-1:m-tau-L+1); % X1 is a Window of X_clean
                %V1 = V1(:); % Vettore colonna
                %PHI_V1 = alpha * PHI_V1 + (1-alpha) * (V1 * V1') / L;
                PHI_V1 = xcorr(V1, V1','biased');
                % Noise correlation vector (1xL)
                gamma_tmp = PHI_V1(1, tau+1 : (tau+L)) / (PHI_V1(1,1)+ 1e-12); % vettore (1x(tau+L-1)) % -1 perché deve essere lungo quanto L!!
                gamma_v = gamma_tmp(:); % ci assicura che sia un vettore colonna
                
                gamma_v_mat(:, k, m) = abs(gamma_v);
            end
        end
        gamma_v_all(tau_idx, noise_idx) = mean(gamma_v_mat, 'all', 'omitnan');
    end
end

%% 3. MEDIA FINALE SU TUTTI I SEGNALI (VETTORI [26 x 1])
gamma_x_mean = mean(gamma_x_all, 2); % Vettore 26x1 per tau = 0..25
gamma_v_mean = mean(gamma_v_all, 2); % Vettore 26x1 per tau = 0..25

end