%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rappresentazione Fig.7 pag.14 paper

% Lettura delle Metriche_%d.mat oppurtunamente calcolate tramite "main3_parfor.m". 
% Lettura delle metriche più la formazione di matricione grosse vsd_all in cui inserire tutti i vsd etc. 
% Ricordiamo che metriche vsd, zita etc... sono 10x1. In cui le righe rappresentano il ritardo tau
% Inoltre ci sono 30 Metriche_%d.mat, perché:
% - 1:10 -- noisy speech con home alarm 
% - 11:20 -- noisy speech con hospital beep1
% - 21:30 -- noisy speech con hospital beep2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; clc;

tau = 1:10; % x-axis

vsd_f_all = zeros(10,30); 
vsd_b_all = zeros(10,30);
vsd_c_all = zeros(10,30);
vsd_np_all = zeros(10,30);
zita_nrH_f_all = zeros(10,30);
zita_nrH_b_all = zeros(10,30);
zita_nrH_c_all = zeros(10,30);
zita_nrH_np_all = zeros(10,30);
PESQ_f_all = zeros(10,30);
PESQ_b_all = zeros(10,30);
PESQ_c_all = zeros(10,30);
PESQ_np_all = zeros(10,30);
STOI_f_all = zeros(10,30);
STOI_b_all = zeros(10,30);
STOI_c_all = zeros(10,30);
STOI_np_all = zeros(10,30);


for speech_index = 1:30 
    metrics_file = sprintf('Metriche_%d.mat',speech_index);
    path_metrics = fullfile('results',metrics_file);
    data = load(path_metrics); %importazione di tutte le metriche
    
    % Memorizzazione di tutte le metriche
    vsd_f_all(:, speech_index) = data.vsd_f;
    vsd_b_all(:, speech_index) = data.vsd_b;
    vsd_c_all(:, speech_index) = data.vsd_c;
    vsd_np_all(:, speech_index) = data.vsd_np;
    zita_nrH_f_all(:, speech_index) = data.zita_nrH_f;
    zita_nrH_b_all(:, speech_index) = data.zita_nrH_b;
    zita_nrH_c_all(:, speech_index) = data.zita_nrH_c;
    zita_nrH_np_all(:, speech_index) = data.zita_nrH_np;
    PESQ_f_all(:, speech_index) = data.PESQ_f;
    PESQ_b_all(:, speech_index) = data.PESQ_b;
    PESQ_c_all(:, speech_index) = data.PESQ_c;
    PESQ_np_all(:, speech_index) = data.PESQ_np;
    STOI_f_all(:, speech_index) = data.STOI_f;
    STOI_b_all(:, speech_index) = data.STOI_b;
    STOI_c_all(:, speech_index) = data.STOI_c;
    STOI_np_all(:, speech_index) = data.STOI_np;
end


%% Calcolo delle metriche medie sui 30 speech (vettori 10x1)
% La media viene fatta lungo la dimensione 2 (le colonne, ossia secondo lo speech_index)
% - 1:10 prendiamo solo i segnali con rumore home alarm --> media salvata in posizione 1 delle martici ..._mean(:,1)
% - 11:20 prendiamo solo i segnali con rumore beep1 ------> media salvata in posizione 2 delle martici ..._mean(:,2)
% - 21:30 prendiamo solo i segnali con rumore beep2 ------> media salvata in posizione 3 delle martici ..._mean(:,3)

%Inizializzazione matrici in cui salveremo la media
vsd_f_mean = zeros(10,3); zita_nrH_f_mean = zeros(10,3); PESQ_f_mean = zeros(10,3); STOI_f_mean = zeros(10,3);
vsd_b_mean = zeros(10,3); zita_nrH_b_mean = zeros(10,3); PESQ_b_mean = zeros(10,3); STOI_b_mean = zeros(10,3);
vsd_c_mean = zeros(10,3); zita_nrH_c_mean = zeros(10,3); PESQ_c_mean = zeros(10,3); STOI_c_mean = zeros(10,3);
vsd_np_mean = zeros(10,3); zita_nrH_np_mean = zeros(10,3); PESQ_np_mean = zeros(10,3); STOI_np_mean = zeros(10,3);

v = 1; %Indice per ciclo for

for i=1:3 % 1 --> speech da 1-10, 

    vsd_f_mean(:,i) = mean(vsd_f_all(:, v:10*i), 2); 
    vsd_b_mean(:,i) = mean(vsd_b_all(:, v:10*i), 2);
    vsd_c_mean(:,i) = mean(vsd_c_all(:, v:10*i), 2);
    vsd_np_mean(:,i) = mean(vsd_np_all(:, v:10*i), 2);

    zita_nrH_f_mean(:,i) = mean(zita_nrH_f_all(:, v:10*i), 2);
    zita_nrH_b_mean(:,i) = mean(zita_nrH_b_all(:, v:10*i), 2);
    zita_nrH_c_mean(:,i) = mean(zita_nrH_c_all(:, v:10*i), 2);
    %zita_nrH_np_mean(:,i) = mean(zita_nrH_np_all(:, v:10*i), 2, 'omitnan');
    zita_nrH_np_mean(:,i) = ones(10,1); % per rendere unprocessed a 0 dB

    PESQ_f_mean(:,i) = mean(PESQ_f_all(:, v:10*i), 2);
    PESQ_b_mean(:,i) = mean(PESQ_b_all(:, v:10*i), 2);
    PESQ_c_mean(:,i) = mean(PESQ_c_all(:, v:10*i), 2);
    PESQ_np_mean(:,i) = mean(PESQ_np_all(:, v:10*i), 2);

    STOI_f_mean(:,i) = mean(STOI_f_all(:, v:10*i), 2);
    STOI_b_mean(:,i) = mean(STOI_b_all(:, v:10*i), 2);
    STOI_c_mean(:,i) = mean(STOI_c_all(:, v:10*i), 2);
    STOI_np_mean(:,i) = mean(STOI_np_all(:, v:10*i), 2);

    v = v + 10; %scorriamo di 10 speech
end


fprintf('Tutte le metriche sono state caricate e mediate con successo!\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grafici %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (pag.14) - Realizzazione 15 grafici con vsd(tau), zita(tau), PESQ(tau), STOI(tau), spettrogrammi e_f(tau),e_b(tau),e_c(tau) con L=3
%     con  tau = [1:10]

%% Plot indici vsd

%Normalizzazione vsd
max_vsd_f = max(vsd_f_mean, [] ,'all');
max_vsd_b = max(vsd_b_mean, [] ,'all');
max_vsd_c = max(vsd_c_mean, [] ,'all');
max_vsd_np = max(vsd_np_mean, [] ,'all');
max_vsd = max([max_vsd_f, max_vsd_b, max_vsd_c, max_vsd_np]);

vsd_f_norm = vsd_f_mean/max_vsd;
vsd_b_norm = vsd_b_mean/max_vsd;
vsd_c_norm = vsd_c_mean/max_vsd;
vsd_np_norm = vsd_np_mean/max_vsd;

% Indice vsd per HOME ALARM
figure;
subplot(1,3,1); %posizionamento dei grafici
plot(tau,10*log10(vsd_f_norm(:,1)),'r-o',tau,10*log10(vsd_b_norm(:,1)),'y-s',tau,10*log10(vsd_c_norm(:,1)),'b-d',tau,10*log10(vsd_np_norm(:,1)),'g','LineWidth', 1.2); 
title("Home Alarm Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("v_{sd}[dB]",'FontWeight','bold');

% Indice vsd per BEEP1
subplot(1,3,2); %posizionamento dei grafici
plot(tau,10*log10(vsd_f_norm(:,2)),'r-o',tau,10*log10(vsd_b_norm(:,2)),'y-s',tau,10*log10(vsd_c_norm(:,2)),'b-d',tau,10*log10(vsd_np_norm(:,2)),'g','LineWidth', 1.2);
title("Beep_{1} Noise");
grid on;
%% 
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("v_{sd}[dB]",'FontWeight','bold');

% Indice vsd per BEEP2
subplot(1,3,3); %posizionamento dei grafici
plot(tau,10*log10(vsd_f_norm(:,3)),'r-o',tau,10*log10(vsd_b_norm(:,3)),'y-s',tau,10*log10(vsd_c_norm(:,3)),'b-d',tau,10*log10(vsd_np_norm(:,3)),'g','LineWidth', 1.2);
title("Beep_{2} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("v_{sd}[dB]",'FontWeight','bold');


%% Plot indici zita_nrH (in decibel)

% Indici zita_nrH per HOME ALARM
figure;
subplot(1,3,1); %posizionamento dei grafici
plot(tau,10*log10(zita_nrH_f_mean(:,1)),'r-o',tau,10*log10(zita_nrH_b_mean(:,1)),'y-s',tau,10*log10(zita_nrH_c_mean(:,1)),'b-d',tau,10*log10(zita_nrH_np_mean(:,1)),'g','LineWidth', 1.2);
title("Home Alarm Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 12]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("\xi_{nr} [dB]",'FontWeight','bold');

% Indici zita_nrH per BEEP1
subplot(1,3,2); %posizionamento dei grafici
plot(tau,10*log10(zita_nrH_f_mean(:,2)),'r-o',tau,10*log10(zita_nrH_b_mean(:,2)),'y-s',tau,10*log10(zita_nrH_c_mean(:,2)),'b-d',tau,10*log10(zita_nrH_np_mean(:,2)),'g','LineWidth', 1.2);
title("Beep_{1} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 12]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("\xi_{nr} [dB]",'FontWeight','bold');

% Indici zita_nrH per BEEP2
subplot(1,3,3); %posizionamento dei grafici
plot(tau,10*log10(zita_nrH_f_mean(:,3)),'r-o',tau,10*log10(zita_nrH_b_mean(:,3)),'y-s',tau,10*log10(zita_nrH_c_mean(:,3)),'b-d',tau,10*log10(zita_nrH_np_mean(:,3)),'g','LineWidth', 1.2);
title("Beep_{2} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 12]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("\xi_{nr} [dB]",'FontWeight','bold');

%% Plot indici PESQ

% Plot per HOME ALARM
figure;
subplot(1,3,1); %posizionamento dei grafici
plot(tau,PESQ_f_mean(:,1),'r-o',tau,PESQ_b_mean(:,1),'y-s',tau,PESQ_c_mean(:,1),'b-d',tau,PESQ_np_mean(:,1),'g','LineWidth', 1.2);
title("Home Alarm Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([1.2 1.9]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("PESQ",'FontWeight','bold');

% Plot per BEEP1
subplot(1,3,2); %posizionamento dei grafici
plot(tau,PESQ_f_mean(:,2),'r-o',tau,PESQ_b_mean(:,2),'y-s',tau,PESQ_c_mean(:,2),'b-d',tau,PESQ_np_mean(:,2),'g','LineWidth', 1.2);
title("Beep_{1} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([1.2 1.9]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("PESQ",'FontWeight','bold');

% Plot per BEEP2
subplot(1,3,3); %posizionamento dei grafici
plot(tau,PESQ_f_mean(:,3),'r-o',tau,PESQ_b_mean(:,3),'y-s',tau,PESQ_c_mean(:,3),'b-d',tau,PESQ_np_mean(:,3),'g','LineWidth', 1.2);
title("Beep_{2} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([1.2 1.9]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("PESQ",'FontWeight','bold');

%% Plot indici STOI

figure;

% Plot per HOME ALARM
subplot(1,3,1); %posizionamento dei grafici
plot(tau,STOI_f_mean(:,1),'r-o',tau,STOI_b_mean(:,1),'y-s',tau,STOI_c_mean(:,1),'b-d',tau,STOI_np_mean(:,1),'g','LineWidth', 1.2);
title("Home Alarm Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0.86 0.98]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("STOI",'FontWeight','bold');

% Plot per BEEP1
subplot(1,3,2); %posizionamento dei grafici
plot(tau,STOI_f_mean(:,2),'r-o',tau,STOI_b_mean(:,2),'y-s',tau,STOI_c_mean(:,2),'b-d',tau,STOI_np_mean(:,2),'g','LineWidth', 1.2);
title("Beep_{1} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0.86 0.98]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("STOI",'FontWeight','bold');

% Plot per HOME ALARM
subplot(1,3,3); %posizionamento dei grafici
plot(tau,STOI_f_mean(:,3),'r-o',tau,STOI_b_mean(:,3),'y-s',tau,STOI_c_mean(:,3),'b-d',tau,STOI_np_mean(:,3),'g','LineWidth', 1.2);
title("Beep_{2} Noise");
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0.86 0.98]); %per limitare asse y
legend('F trad','B trad','C trad','unprocessed',Location ='best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("STOI",'FontWeight','bold');

%% Spettrogrammi RUMORI

Fs = 16000;

% --- Parametri STFT ---
win_length = 512;
win = hamming(win_length); 
overlap = round(win_length * 0.75); 
nfft = 512;

% Prendiamo i rumori (dominio freq)
[V_home_alarm, V_beep1, V_beep2] = my_all_noises(Fs);
% E li riportiamo nel dominio del tempo
home_alarm = istft(V_home_alarm, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
beep1 = istft(V_beep1, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
beep2 = istft(V_beep2, Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);

figure;
subplot(1,3,1); %posizionamento dei grafici
spectrogram(home_alarm, win, overlap, nfft, Fs, 'yaxis');
title('Home Alarm Noise');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontWeight', 'bold');

subplot(1,3,2); %posizionamento dei grafici
spectrogram(beep1, win, overlap, nfft, Fs, 'yaxis');
title('Hospital Beep 1 Noise');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontWeight', 'bold');

subplot(1,3,3); %posizionamento dei grafici
spectrogram(beep2, win, overlap, nfft, Fs, 'yaxis');
title('Hospital Beep 2 Noise');
xlabel('Time [s]', 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontWeight', 'bold');