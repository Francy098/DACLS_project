%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcolo metriche per grafici Fig.9 pag.16
%
% Qui ci occupiamo di  determinarci lo NRR (Noise Reduction Ratio) per
% frame ad 1kHz filtrando opportunatamente un speech di voce femminile che
% dura 3,4s corrotto con il rumore Beep1 ad un SNR di 10dB. Ciò che
% vogliamo rappresentare sono le seguenti condizioni:
% - NRR per L = 3, tau = 3, step_size = 0.5, idx_threshold = -25 dB
% - NRR per L = 3, tau = 3, step_size = 0.5, idx_threshold = -40 dB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Inizializzazione parametri

clear all; clc; 

L = 3; % tappi filtro
tau = 3; %Campioni di ritardo

Fs = 16000; % 16kHz
%tau=10; campioni di ritardo

% Parametri stft
win_length = 512;
win = hamming(win_length); 
% Overlap
overlap = round(win_length * 0.75); %per renderlo al 75%
nfft = 512;


%% Lettura Speech rumorosi e PULITI  

[y_20k, fs1] = audioread('sounds/noisy_female_Optimal_SNR10dB.wav'); % speech rumoroso
[w_20k, fs2] = audioread('our_TIMIT_database/female_Optimal.wav'); % speech pulito

% Ricampionamento a 16 kHz
y_noisy = resample(y_20k, Fs, fs1);
y_noisy = y_noisy(:); % Assicura che sia un vettore colonna [N x 1]
x_clean = resample(w_20k, Fs, fs2);  
x_clean = x_clean(:); % Assicura che sia un vettore colonna [N x 1]

% STFT
M_target = 449;
X_clean_temp = stft(x_clean, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap); 
Y1_temp      = stft(y_noisy, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap); 
% Rideminionamento per precauzione
X_clean = X_clean_temp(:,1:M_target);
Y1 = Y1_temp(:,1:M_target);

[K, M] = size(Y1); % K = frequenze, M = frame totali

% Indice del bin di frequenza corrispondente a 1 kHz (1000 Hz)
% Freq_bin = round(f_target / (Fs / nfft)) + 1
%f_bin_1kHz = round(1000 / (Fs / nfft)) + 1;

% Preallocazione matrici NRR
NRR_f = zeros(2, M);
NRR_b = zeros(2, M);
NRR_c = zeros(2, M);

%% Lettura rumori 
% Ci procuriamo solo il rumore Beep1
[~, V_noise, ~,] = my_three_noises(Fs);
  
    
%% Preallocazione matrici
Hw_f = zeros(L, K, M);
Hw_b = zeros(L, K, M);
Hw_c = zeros(L, K, M);


%% INIZIO FILTRAGGIO (idx_threshold = -25dB)
fprintf('Start filtering!\n');

    
%Adaptive Filtering 
[ ~,~,Hw_f, ~,~,Hw_b, ~,~,Hw_c, Lf_bank,Lb_bank,Lc_bank] = Adaptive_Filtering(tau,L,Y1,V_noise);
        

% Performance measures
[~,~, nrr_f_full] = my_performance_measures(tau,L,Hw_f,V_noise,X_clean,Lf_bank); %Forward
[~,~, nrr_b_full] = my_performance_measures(tau,L,Hw_b,V_noise,X_clean,Lb_bank); %Backward
[~,~, nrr_c_full] = my_performance_measures(tau,L,Hw_c,V_noise,X_clean,Lc_bank); %Combined Approch

% nrr _%_full è una martrice che contiene tutti i bin delle frequenza. Noi
% ci prendiamo solo quella che ci interessa

NRR_f(1,:) = nrr_f_full;%(f_bin_1kHz, :);
NRR_b(1,:) = nrr_b_full;%(f_bin_1kHz, :);
NRR_c(1,:) = nrr_c_full;%(f_bin_1kHz, :);


%% INIZIO FILTRAGGIO (idx_threshold = -40 dB)
%fprintf('Start filtering!\n');
    
    
%Adaptive Filtering 
[ ~,~,Hw_f, ~,~,Hw_b, ~,~,Hw_c, Lf_bank,Lb_bank,Lc_bank] = Adaptive_Filtering_40idx(tau,L,Y1,V_noise);
        

% Performance measures
[~,~, nrr_f_full] = my_performance_measures(tau,L,Hw_f,V_noise,X_clean,Lf_bank); %Forward
[~,~, nrr_b_full] = my_performance_measures(tau,L,Hw_b,V_noise,X_clean,Lb_bank); %Backward
[~,~, nrr_c_full] = my_performance_measures(tau,L,Hw_c,V_noise,X_clean,Lc_bank); %Combined Approch

NRR_f(2,:) = nrr_f_full;%(f_bin_1kHz, :);
NRR_b(2,:) = nrr_b_full;%(f_bin_1kHz, :);
NRR_c(2,:) = nrr_c_full;%(f_bin_1kHz, :);

fprintf('End filtering!\n');

%% Salvataggio metriche
save('results/Metriche_Adaptive_NRR.mat', 'NRR_f', 'NRR_b', 'NRR_c');  
fprintf('Saved in Metriche_Adaptive_NRR.mat\n');

 
%% Grafici

load('results/Metriche_Adaptive_NRR.mat'); % caricamento metriche

frames = 1:M; 
interval = 1:M; 

figure('Name', 'Fig.9 - NRR Comparison', 'NumberTitle', 'off');

% Aggiustamento offset
NRR_f(1,:) = NRR_f(1,:) +25; % +25dB
NRR_b(1,:) = NRR_b(1,:) +25; % +25dB
NRR_c(1,:) = NRR_c(1,:) +25; % +25dB

% Subplot (a) per idx_threshold = -25 dB
subplot(1, 2, 1);
plot(frames(interval), NRR_f(1, interval), 'r-o', 'LineWidth', 1.2); hold on;
plot(frames(interval), NRR_b(1, interval), 'y-s', 'LineWidth', 1.2);
plot(frames(interval), NRR_c(1, interval), 'b-d', 'LineWidth', 1.2);
%yline(spectral_level_1kHz_dB, 'm--', 'LineWidth', 1.2, 'DisplayName', 'Spectral noise floor (+25 dB)');
ylim([-40 40]);
xlim([170 261]);
hold off; grid on;
title('(a) thres = -25dB', 'FontWeight', 'bold');
ylabel('NRR [dB]', 'FontWeight', 'bold');
xlabel('Frame', 'FontWeight', 'bold');
legend('Filter F', 'Filter B', 'Filter C', 'Location', 'best');

% Aggiustamento offset
NRR_f(2,:) = NRR_f(2,:) +20; % +20dB
NRR_b(2,:) = NRR_b(2,:) +20; % +20dB
NRR_c(2,:) = NRR_c(2,:) +20; % +20dB

% Subplot (b) per idx_threshold = -40 dB
subplot(1, 2, 2);
plot(frames(interval), NRR_f(2, interval), 'r-o', 'LineWidth', 1.2); hold on;
plot(frames(interval), NRR_b(2, interval), 'y-s', 'LineWidth', 1.2);
plot(frames(interval), NRR_c(2, interval), 'b-d', 'LineWidth', 1.2);
%yline(spectral_level_1kHz_dB, 'm--', 'LineWidth', 1.2, 'DisplayName', 'Spectral noise floor (+25 dB)');
ylim([-40 40]);
xlim([170 261]);
hold off; grid on;
title('(b) thres = -40dB', 'FontWeight', 'bold');
ylabel('NRR [dB]', 'FontWeight', 'bold');
xlabel('Frame', 'FontWeight', 'bold');
legend('Filter F', 'Filter B', 'Filter C', 'Location', 'best');