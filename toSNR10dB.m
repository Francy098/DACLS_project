%% Caricamento speech e rumore

Fs=16000;

for i = 1:5

%path_speech = sprintf('our_TIMIT_database/male_%d.wav',i);
path_speech = sprintf('our_TIMIT_database/female_%d.wav',i);
[x_20k, fs] = audioread(path_speech); 
x = resample(x_20k, Fs, fs);

% Caricamento rumore
%[v_20k, fs] = audioread("rumori\House_alarm\home_alarm1.flac");
[v_20k, fs] = audioread("rumori\Heart_monitor_beeping\medical_noise2.wav");
v = resample(v_20k, Fs, fs);

% Se il rumore è stereo, lo convertiamo in mono
if size(v, 2) > 1
    v = mean(v, 2);
end

% Uniformiamo la lunghezza del rumore a quella dello speech (N campioni)
N = length(x);
if length(v) >= N
    v = v(1:N); % Tagliamo il rumore in eccesso
else
    v = [v; zeros(N - length(v), 1)]; % Applica zero-padding se il rumore è più corto
end

%% Calcolo del guadagno per SNR = 10 dB
% Calcolo delle potenze dei due segnali originari
P_x = mean(x.^2); % Potenza dello speech
P_v = mean(v.^2); % Potenza del rumore non scalato

% Calcolo il fattore di scala g per il rumore
SNR_dB = 10;
g = sqrt(P_x / (P_v * 10^(SNR_dB / 10)));

% Somma lo speech e il rumore scalato
y1 = x + g * v;

%path_noisy_speech = sprintf('sounds/home_alarm/noisy_female_%d_SNR10dB.wav',i);
%path_noisy_speech = sprintf('sounds/medical_beep1/noisy_female_%d_SNR10dB.wav',i);
path_noisy_speech = sprintf('sounds/medical_beep2/noisy_female_%d_SNR10dB.wav',i);
audiowrite(path_noisy_speech, y1, Fs);

end