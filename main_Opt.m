%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcolo e rappresentazione metriche del FILTRO OTTIMO (pag.11-12)

% In questo script ci focalizziamo sull'implementazione del filtro OTTIMO e
% alla rappresentazione dei grafici del paper presenti nell'articolo a pag.11-12

% Il rumore è composto dalla generazione di più armoniche modulute tramite
% un treno trapezoidale nel dominio del tempo con il rumore bianco gaussiano
% noise detector -25dB

%Fig.3 Spettrogramma di un segnale lungo 3.4s famale speech corrotto da un
%rumore sintetico non stazionario a SNR = 10dB

%Fig.4 (vsd, zita_nr, PESQ e STOI) Female noisy speech di 3.4s, filtrato con LS ALE filter sempre a L=3, variando tau [1-10]

%Fig.5 (vsd, zita_nr) Female noisy speech di 3.4s, filtrato con CLP (combined approach) LS


%% Calcolo e rappresentazione metriche del FILTRO OTTIMO (pag.11-12)
% In questo script ci focalizziamo sull'implementazione del filtro OTTIMO e
% alla rappresentazione dei grafici del paper presenti nell'articolo a pag.11-12.

clear; clc; close all;

%% Inizializzazione parametri
Fs = 16000; % 16 kHz

% Parametri STFT
win_length = 512;
win = hamming(win_length); 
overlap = round(win_length * 0.75); % Overlap al 75%
nfft = 512;

% Configurazione PESQ
binary = 'pesq2.exe';
pathaudio = 'sounds';
if ~exist(pathaudio, 'dir')
    mkdir(pathaudio);
end

reference_f = 'filtered_f_opt.wav';
reference_b = 'filtered_b_opt.wav';
reference_c = 'filtered_c_opt.wav';

% Preallocazione metriche scalari (funzione di tau = 1:10)
vsd_f = zeros(10,1); zita_nrH_f = zeros(10,1); PESQ_f = zeros(10,1); STOI_f = zeros(10,1);
vsd_b = zeros(10,1); zita_nrH_b = zeros(10,1); PESQ_b = zeros(10,1); STOI_b = zeros(10,1);
vsd_c = zeros(10,1); zita_nrH_c = zeros(10,1); PESQ_c = zeros(10,1); STOI_c = zeros(10,1);

%% Lettura Speech rumorosi e puliti
[y_20k, fs1] = audioread('sounds/noisy_female_Optimal_SNR10dB.wav');
[w_20k, fs2] = audioread('our_TIMIT_database/female_Optimal.wav');

% Ricampionamento a 16 kHz
y_noisy = resample(y_20k, Fs, fs1);
y_noisy = y_noisy(:); % Assicura che sia un vettore colonna [N x 1]
x_clean = resample(w_20k, Fs, fs2);  
x_clean = x_clean(:); % Assicura che sia un vettore colonna [N x 1]

% Caricamento rumore di riferimento
[~, ~, V_noise] = my_three_noises(Fs); 

% STFT
X_clean = stft(x_clean, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap); 
Y1      = stft(y_noisy, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap); 

[K, M] = size(Y1); % K = frequenze, M = frame totali

% Preallocazione matrici NRR
NRR_f = zeros(2, M);
NRR_b = zeros(2, M);
NRR_c = zeros(2, M);

% Generazione e allineamento per ISTFT e STOI
y1 = istft(Y1, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
len_istft = length(y1);

audiowrite('sounds/degraded_Opt.wav', y_noisy(1:len_istft), Fs); 
degraded_file = 'degraded_Opt.wav'; %for PESQ
reference_stoi = x_clean(1:len_istft);

% --- Generazione del file reference pulito temporaneo nella cartella 'sounds' ---
clean_file = sprintf('clean_speech_Opt.wav');
clean_path = fullfile('sounds', clean_file);
audiowrite(clean_path, reference_stoi, Fs);

% Preallocazione segnali tempo
e_l = zeros(len_istft, 10);
e_f = zeros(len_istft, 10);
e_b = zeros(len_istft, 10);
e_c = zeros(len_istft, 10);

L = 3; % Taps del filtro

% Creazione filtro Hw_np --> per il rappresentazione metriche del non processato
% Siccome E = Y1, basta implementare un Hw_np pieno di zeri
Hw_np = zeros(L, K, M);

% Inizializziamo le sue relative metriche che ci interessano graficare
vsd_np     = zeros(10,1); 
zita_nrH_np = zeros(10,1);
PESQ_np = zeros(10,1); 
STOI_np = zeros(10,1);

%% Optimal Filtering Loop
fprintf('Start filtering!\n');

for tau = 1:10
    fprintf('tau = %d...\n', tau);
    
    % Esecuzione Filtraggio Ottimo
    [E_f, Z_f, Hw_f, E_b, Z_b, Hw_b, E_c, Z_c, Hw_c] = Optimal_Filtering(tau, L, Y1, V_noise);
    
    % Ricostruzione dominio del tempo
    e_f(:, tau) = istft(E_f, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
    e_b(:, tau) = istft(E_b, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
    e_c(:, tau) = istft(E_c, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
    
    % Scrittura audio per PESQ
    audiowrite('sounds/filtered_f_opt.wav', e_f(:, tau), Fs);
    audiowrite('sounds/filtered_b_opt.wav', e_b(:, tau), Fs);
    audiowrite('sounds/filtered_c_opt.wav', e_c(:, tau), Fs); 
    
    %% Performance Measures
    if tau == 1
        [zita_nrH_f(tau), vsd_f(tau), NRR_f(1, :)] = Optimal_Performance_Measures(tau, L, Hw_f, V_noise, X_clean);
        [zita_nrH_b(tau), vsd_b(tau), NRR_b(1, :)] = Optimal_Performance_Measures(tau, L, Hw_b, V_noise, X_clean);
        [zita_nrH_c(tau), vsd_c(tau), NRR_c(1, :)] = Optimal_Performance_Measures(tau, L, Hw_c, V_noise, X_clean);
        [zita_nrH_np(tau), vsd_np(tau), ~] = Optimal_Performance_Measures(tau, L, Hw_np, V_noise, X_clean);
    elseif tau == 3
        [zita_nrH_f(tau), vsd_f(tau), NRR_f(2, :)] = Optimal_Performance_Measures(tau, L, Hw_f, V_noise, X_clean);
        [zita_nrH_b(tau), vsd_b(tau), NRR_b(2, :)] = Optimal_Performance_Measures(tau, L, Hw_b, V_noise, X_clean);
        [zita_nrH_c(tau), vsd_c(tau), NRR_c(2, :)] = Optimal_Performance_Measures(tau, L, Hw_c, V_noise, X_clean);
        [zita_nrH_np(tau), vsd_np(tau), ~] = Optimal_Performance_Measures(tau, L, Hw_np, V_noise, X_clean);
    else
        [zita_nrH_f(tau), vsd_f(tau), ~] = Optimal_Performance_Measures(tau, L, Hw_f, V_noise, X_clean);
        [zita_nrH_b(tau), vsd_b(tau), ~] = Optimal_Performance_Measures(tau, L, Hw_b, V_noise, X_clean);
        [zita_nrH_c(tau), vsd_c(tau), ~] = Optimal_Performance_Measures(tau, L, Hw_c, V_noise, X_clean);
        [zita_nrH_np(tau), vsd_np(tau), ~] = Optimal_Performance_Measures(tau, L, Hw_np, V_noise, X_clean);
    end
    
    
    % Valutazione PESQ
    PESQ_f(tau) = pesq2_mtlb(clean_file, reference_f, Fs, 'wb', binary, pathaudio); %forward
    PESQ_b(tau) = pesq2_mtlb(clean_file, reference_b, Fs, 'wb', binary, pathaudio); %backward
    PESQ_c(tau) = pesq2_mtlb(clean_file, reference_c, Fs, 'wb', binary, pathaudio); %combined
    PESQ_np(tau) = pesq2_mtlb(clean_file, degraded_file, Fs, 'wb', binary, pathaudio); %unprocessed
    
    % Valutazione STOI
    STOI_f(tau) = stoi(e_f(:, tau), reference_stoi, Fs); %forward
    STOI_b(tau) = stoi(e_b(:, tau), reference_stoi, Fs); %backward
    STOI_c(tau) = stoi(e_c(:, tau), reference_stoi, Fs); %cambined
    STOI_np(tau) = stoi(y_noisy, reference_stoi, Fs); %unprocessed

end

fprintf('End filtering!\n');


%% Salvataggio risultati

save('results\Metriche_Opt.mat', ...
    'vsd_f', 'zita_nrH_f', 'PESQ_f', 'STOI_f', ...
    'vsd_b', 'zita_nrH_b', 'PESQ_b', 'STOI_b', ...
    'vsd_c', 'zita_nrH_c', 'PESQ_c', 'STOI_c', ...
    'vsd_np', 'zita_nrH_np', 'PESQ_np', 'STOI_np', ...
    'NRR_f', 'NRR_b', 'NRR_c');
fprintf('Metrics saved!\n');

%% Calculate gamma_x e gamma_v for all speeches and noises
%
% In questo caso abbiamo utilizzato l'intero database condiderato nel progetto:
% - 20 speech (10 female voice + 10 male voice)
% - 26 nonstationary harmonic noise signals (heart monitor beeping, train door beeping, house alarm, smoke alarm, and rail road crossing bells)
%
% I gamma vengono calcolati per ogni speech (gamma_x) e rumori (gamma_v). 
% Dopo di che facciamo una media di tutti i gamma_x e gamma_v per trovarne un unico rappresentativo
% del gruppo.

fprintf('Start gamma calculation...\n');
[gamma_x, gamma_v] = new_gamma_calculation(K, M, L, Fs);
fprintf('End gamma calculation!\n');
% Salvataggio gamma_x e gamma_v
save('results\gamma_calculation.mat','gamma_x', 'gamma_v')
fprintf('Metrics saved!\n');

% DA FARE:
% - creazione di una my_all_database, che restituisce un due matrici: una
% con tutti gli speech ed un'altra con tutti i rumori
% - adattare gamma_calculation.m per più segnali



%% Grafici: Spettrogramma e Coefficienti di Correlazione

tau_x_axis = 0:25;

load('results\gamma_calculation.mat'); % caricamento gamma

% Normalizzazione
%max_gamma = max(max(gamma_x),max(gamma_v));
%gamma_x = gamma_x ./ max(gamma_x);
%gamma_v = gamma_v ./ max(gamma_v);

figure('Name', 'Spectrogram and Correlation', 'NumberTitle', 'off');
subplot(1, 2, 1);
spectrogram(y_noisy, win, overlap, nfft, Fs, 'yaxis');
title('Spectrogram of Noisy Speech');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontWeight', 'bold');

subplot(1, 2, 2);
plot(tau_x_axis, abs(gamma_x), 'b-o', 'LineWidth', 1.2); hold on;
plot(tau_x_axis, abs(gamma_v), 'r-s', 'LineWidth', 1.2); hold off;
grid on;
%ylim([0 1]);
legend('Clean speech signals', 'Non-stationary harmonic noise signals', 'Location', 'best');
xlabel('\tau', 'FontWeight', 'bold');
ylabel('|\gamma|', 'FontWeight', 'bold');
title('Correlation Coefficients Magnitude');


%% Grafici: Indici di Qualità

load('results\Metriche_Opt.mat'); % caricamento metriche

tau_x_axis = 1:10;

figure('Name', 'Quality Metrics Comparison', 'NumberTitle', 'off');

%Normalizzazione vsd
max_vsd_f = max(vsd_f, [] ,'all');
max_vsd_b = max(vsd_b, [] ,'all');
max_vsd_c = max(vsd_c, [] ,'all');
max_vsd_np = max(vsd_np, [] ,'all');
max_vsd = max([max_vsd_f, max_vsd_b, max_vsd_c, max_vsd_np]);

vsd_f_norm = vsd_f/max_vsd;
vsd_b_norm = vsd_b/max_vsd;
vsd_c_norm = vsd_c/max_vsd;
vsd_np_norm = vsd_np/max_vsd;

% Spectral Distortion Index (v_sd)
subplot(2, 2, 1);
plot(tau_x_axis, vsd_f_norm, 'r-o', tau_x_axis, vsd_b_norm, 'y-s', tau_x_axis, vsd_c_norm, 'b-d', tau_x_axis,  vsd_np_norm, 'g', 'LineWidth', 1.2);
%ylim([0 1]);
grid on;
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel('\tau', 'FontWeight', 'bold');
ylabel('v_{sd}', 'FontWeight', 'bold');
title('Spectral Distortion Index (v_{sd})');

% Noise Reduction Factor (\xi_{nr})
subplot(2, 2, 2);
plot(tau_x_axis, 10*log10(zita_nrH_f), 'r-o', tau_x_axis, 10*log10(zita_nrH_b), 'y-s', tau_x_axis, 10*log10(zita_nrH_c), 'b-d', tau_x_axis,  10*log10(zita_nrH_np), 'g', 'LineWidth', 1.2);
grid on;
legend('F trad','B trad','C trad','unprocessed',Location ='best');
xlabel('\tau', 'FontWeight', 'bold');
ylabel('\xi_{nr} [dB]', 'FontWeight', 'bold');
title('Noise Reduction Factor (\xi_{nr})');

% PESQ Score
subplot(2, 2, 3);
plot(tau_x_axis, PESQ_f, 'r-o', tau_x_axis, PESQ_b, 'y-s', tau_x_axis, PESQ_c, 'b-d', tau_x_axis, PESQ_np, 'g', 'LineWidth', 1.2);
grid on;
legend('F trad','B trad','C trad','unprocessed',Location ='best');
xlabel('\tau', 'FontWeight', 'bold');
ylabel('PESQ', 'FontWeight', 'bold');
title('PESQ Score');

% STOI Score
subplot(2, 2, 4);
plot(tau_x_axis, STOI_f, 'r-o', tau_x_axis, STOI_b, 'y-s', tau_x_axis, STOI_c, 'b-d', tau_x_axis, STOI_np,'g', 'LineWidth', 1.2);
grid on;
legend('F trad','B trad','C trad','unprocessed',Location ='best');
xlabel('\tau', 'FontWeight', 'bold');
ylabel('STOI', 'FontWeight', 'bold');
title('STOI Score'); 


%% Rappresentazione grafici Fig. 6 (pag. 13 del paper)
frames = 1:M; 
interval = 1:M; 

figure('Name', 'Figura 6 - NRR Comparison', 'NumberTitle', 'off');

% Subplot (a) per tau = 1
subplot(1, 2, 1);
plot(frames(interval), NRR_f(1, interval), 'r-o', 'LineWidth', 1.2); hold on;
plot(frames(interval), NRR_b(1, interval), 'y-s', 'LineWidth', 1.2);
plot(frames(interval), NRR_c(1, interval), 'b-d', 'LineWidth', 1.2);
%yline(spectral_level_1kHz_dB, 'm--', 'LineWidth', 1.2, 'DisplayName', 'Spectral noise floor (+25 dB)');
ylim([-60 60]);
xlim([180 290]);
hold off; grid on;
title('(a) \tau = 1', 'FontWeight', 'bold');
ylabel('NRR [dB]', 'FontWeight', 'bold');
xlabel('Frame', 'FontWeight', 'bold');
legend('Filter F', 'Filter B', 'Filter C', 'Location', 'best');

% Subplot (b) per tau = 3
subplot(1, 2, 2);
plot(frames(interval), NRR_f(2, interval), 'r-o', 'LineWidth', 1.2); hold on;
plot(frames(interval), NRR_b(2, interval), 'y-s', 'LineWidth', 1.2);
plot(frames(interval), NRR_c(2, interval), 'b-d', 'LineWidth', 1.2);
%yline(spectral_level_1kHz_dB, 'm--', 'LineWidth', 1.2, 'DisplayName', 'Spectral noise floor (+25 dB)');
ylim([-60 60]);
xlim([180 290]);
hold off; grid on;
title('(b) \tau = 3', 'FontWeight', 'bold');
ylabel('NRR [dB]', 'FontWeight', 'bold');
xlabel('Frame', 'FontWeight', 'bold');
legend('Filter F', 'Filter B', 'Filter C', 'Location', 'best');