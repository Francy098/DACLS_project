function [gamma_x_mean, gamma_v_mean] = gamma_calculation(K, M, L, Fs)

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
        cross_x = zeros(K, L);
        gamma_x_mat = NaN(L, K, M); % Preallocazione a NaN per ignorare i frame non calcolati
        
        for m = (tau + L):M
            for k = 1:K
                % --- Stima varianza temporale ---
                phi_X(k, m) = alpha * phi_X(k, m-1) + (1-alpha) * abs(X_clean(k, m))^2;
                
                % --- Vettore gamma_x ---
                gamma_x_vec = zeros(L, 1);
                for r = 0:(L - 1)
                    delay_idx = m - (tau + r) + 1;
                    if delay_idx >= 1
                        cross_x(k, r+1) = alpha * cross_x(k, r+1) + (1-alpha) * X_clean(k, m) * conj(X_clean(k, delay_idx));
                        gamma_x_vec(r+1) = cross_x(k, r+1) / (phi_X(k, m) + 1e-6);
                        
                        % Clipping stabilità
                        if abs(gamma_x_vec(r+1)) > 1
                            gamma_x_vec(r+1) = gamma_x_vec(r+1) / abs(gamma_x_vec(r+1));
                        end
                    end
                end
                
                gamma_x_mat(:, k, m) = abs(gamma_x_vec);
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
        cross_v = zeros(K_noise, L);
        gamma_v_mat = NaN(L, K_noise, M_noise);
        
        for m = (tau + L):M_noise
            for k = 1:K_noise % Cicla fino a K_noise, NON K
                phi_V(k, m) = alpha * phi_V(k, m-1) + (1-alpha) * abs(V_noise(k, m))^2;
                
                % --- Vettore gamma_v ---
                gamma_v_vec = zeros(L, 1);
                for r = 0:(L - 1)
                    delay_idx = m - (tau + r) + 1;
                    
                    % Protezione esplicita sugli indici di memoria
                    if ((delay_idx >= 1) && (delay_idx <= M_noise))
                        cross_v(k, r+1) = alpha * cross_v(k, r+1) + (1-alpha) * V_noise(k, m) * conj(V_noise(k, delay_idx));
                        gamma_v_vec(r+1) = cross_v(k, r+1) / (phi_V(k, m) + 1e-6);
                        
                        if abs(gamma_v_vec(r+1)) > 1
                            gamma_v_vec(r+1) = gamma_v_vec(r+1) / abs(gamma_v_vec(r+1));
                        end
                    end
                end
                
                gamma_v_mat(:, k, m) = abs(gamma_v_vec);
            end
        end
        gamma_v_all(tau_idx, noise_idx) = mean(gamma_v_mat, 'all', 'omitnan');
    end
end

%% 3. MEDIA FINALE SU TUTTI I SEGNALI (VETTORI [26 x 1])
gamma_x_mean = mean(gamma_x_all, 2); % Vettore 26x1 per tau = 0..25
gamma_v_mean = mean(gamma_v_all, 2); % Vettore 26x1 per tau = 0..25

end