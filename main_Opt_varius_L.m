%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMAL FILTERING SOLO CON APPROCCIO COMBINATO (Fig.5 a pag. 12)

% In questo script ci occuperemo di riprodurre i grafici in a più valori di L = {3,5,8,10} su ogni tau = [1-10].
% Fig. 5. (a) Distortion index and (b) noise reduction factor for CLP LS ALE filter applied on the 3.4 s female speech signal corrupted by synthetic nonstationary harmonic
% noise at 10 dB SNR for varying delay t and different filter lengths L = {3,5,8,10}. Solid lines are for traditional performance measures definition, the dashed
% lines are for the OD definition in Section 3. Dashed green is for the unprocessed noisy signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inizializzazione parametri
clear all; clc;

Fs = 16000; % 16kHz
%tau=10; campioni di ritardo

% Parametri stft
win_length = 512;
win = hamming(win_length); 
% Overlap
overlap = round(win_length * 0.75); %per renderlo al 75%
nfft = 512;

% Inizializzazione per il calcolo misura PESQ
% name of executable file for PESQ calculation
binary = 'pesq2.exe';
% specify path to folder with reference and degraded audio files in it
pathaudio = 'sounds';
% specify reference (processed) ... 
reference_c = 'filtered_c_opt.wav';
% ... and degraded audio files
degraded = 'degraded_Opt.wav';
%Per la PESQ measure, ogni file .wav vengono riscritti ad ogni ciclo di tau


%% Lettura Speech rumorosi 

% Campionamento noisy speech
[y_20k,fs1] = audioread('sounds/noisy_female_Optimal_SNR10dB.wav');

% Campionamento segnale pulito (per STOI measure)
[w_20k,fs2] = audioread('our_TIMIT_database/female_Optimal.wav');
                                   % fs=20kHz di default

% Ricampionamento del segnale 16 kHz (Fs)
y_noisy = resample(y_20k, Fs, fs1);
x_clean = resample(w_20k, Fs, fs2);  

% STFT dello speech PULITO corrente, quindi NON RUMOROSO. Ci servirà in my_sperformance_measure
X_clean = stft(x_clean,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap); 

% Lettura V_noise
[~, ~, V_noise] = my_three_noises(Fs); % già in dominio della freq.

%Calcolo STFT di un noisy speech alla volta
%STFT -> otteniamo una Y(k,m). cioè un matrice 2D con k = indice freq e m = indice frame
Y1 = stft(y_noisy,Fs,"Window",win,"FFTLength",nfft); %ogni Fs c'è un valore sul vettore
[K, M] = size(Y1); % K = indice di frequenza, M = numero totale di frame

% Ottieniamo la lunghezza esatta dell'istft per questo specifico segnale
y1 = istft(Y1, Fs, "Window", win, "FFTLength", nfft);
len_istft = length(y1);

% Generazione del file degraded.wav per la PESQ
audiowrite('sounds\degraded_Opt.wav', y_noisy, Fs); 
% Per la misura STOI
reference_stoi = x_clean(1:len_istft); % Taglia alla lunghezza coerente con l'istft


% Preallocazione dinamica (per evitare mismatch di dimensioni tra segnali diversi)
e_c = zeros(len_istft, 10);

% Tappi filtro
L_vec = [3;5;8;10];

% Preallocazione metriche
vsd_c = zeros(10,4); % 4 colonne perché L_vec è lungo 4
zita_nrH_c = zeros(10,4); 
PESQ_c = zeros(10,4); 
STOI_c = zeros(10,4);

%% Optimal Filtering

fprintf('Start filtering!\n');

for i = 1:4 % 4 prché L_vec è lungo 4

    L_curr= L_vec(i);
    fprintf('Elaborazione per L = %d...\n', L_curr);
    
    for tau = 1:10
    
    fprintf('tau = %d...\n',tau);

    % Optimal Filtering 
    [~,~,~, ~,~,~, E_c,Z_c,Hw_c] = Optimal_Filtering(tau,L_curr,Y1,V_noise);

    e_c(:,tau) = istft(E_c,Fs,"Window",win,"FFTLength",nfft);
   
    audiowrite('sounds\filtered_c_opt.wav',e_c(:,tau),Fs); % ci servono i .wav per il calcolo del PESQ

    %% Performance measures
    [zita_nrH_c(tau,i),vsd_c(tau,i), ~] = Optimal_Performance_Measures(tau,L_curr,Hw_c,V_noise,X_clean); %Combined Approch
    
    % PESQ measure
    %audiowrite('sounds\degraded.wav',y_all_noisy_speeches(:,speech_index),Fs); %Creazione degraded.wav corrente
    PESQ_c(tau,i) = pesq2_mtlb(reference_c, degraded, Fs, 'wb', binary, pathaudio ); %full-band

    %STOI measure
    % Specify processed... 
    processed_c = e_c(:,tau);
    % reference_stoi definito sopra
    STOI_c(tau,i) = stoi(processed_c, reference_stoi, Fs);

    end

end

fprintf('End filtering!\n');

%% Salvataggio metriche
    
    save('results/Metriche_Opt_varius_L.mat', ...
        'vsd_c', 'zita_nrH_c', 'PESQ_c', 'STOI_c');

    fprintf('Metrics saved!\n');
                                                           
%% Grafici

tau_x_axis = 1 : 10; % asse x

load('results/Metriche_Opt_varius_L.mat'); % caricamento metriche

% ----------- Inizializzazione metriche per il non processato -------------
% Come visto in precedenza, nel main3_parfor.m e con Grafici.m, il non
% processato presenta un vsd = 0 e con uno zita_dB = 0 dB. Quindi creiamo
% due vettori affinché vengano replicate queste condizioni

vsd_np = zeros(1, 10);
zita_np_dB = zeros(1, 10);

%Normalizzazione vsd
vsd_c_norm = vsd_c/max(vsd_c, [] ,'all');

% % Fig. 5(a) Distortion Index vsd
figure;
subplot(2,2,1); %posizionamento dei grafici
plot(tau_x_axis, vsd_c_norm(:,1),'r-o',tau_x_axis, vsd_c_norm(:,2),'y-s',tau_x_axis, vsd_c_norm(:,3),'b-d',tau_x_axis, vsd_c_norm(:,4),'v-',tau_x_axis,vsd_np,'g','LineWidth', 1.2);
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('L= 3','L= 5','L= 8','L= 10','unprocessed',Location= 'best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("v_{sd}",'FontWeight','bold');
title('Fig. 5(a) Distortion Index (v_{sd})');

% Portiamo tutto in dB
zita_nrH_c = 10*log10(abs(zita_nrH_c));

% Fig. 5(b) Noise Reduction Factor zita_nr (in dB)
subplot(2,2,2); %posizionamento dei grafici
plot(tau_x_axis,zita_nrH_c(:,1),'r-o',tau_x_axis,zita_nrH_c(:,2),'y-s',tau_x_axis,zita_nrH_c(:,3),'b-d',tau_x_axis,zita_nrH_c(:,4),'v-',tau_x_axis,zita_np_dB,'g','LineWidth', 1.2);
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('L= 3','L= 5','L= 8','L= 10','unprocessed',Location= 'best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
title('Fig. 5(b) Noise Reduction Factor (\xi_{nr})');
ylabel("\xi_{nr}[dB]",'FontWeight','bold');