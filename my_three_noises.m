function [V_home_alarm, V_beep1, V_beep2] = my_three_noises(Fs)
N_campioni = 60000; 
all_noises = zeros(N_campioni, 3); 

% --- 1. Home Alarm ---
[v_20k, fs] = audioread("rumori\House_alarm\home_alarm1.flac");
v = resample(v_20k, Fs, fs);
if size(v, 2) > 1, v = mean(v, 2); end
v = v(:); % Forzatura vettore colonna

if length(v) >= N_campioni
    all_noises(:,1) = v(1:N_campioni);
else
    all_noises(:,1) = [v; zeros(N_campioni - length(v), 1)];
end

% --- 2. Medical Beep 1 ---
[v1_20k, fs] = audioread("rumori\Heart_monitor_beeping\medical_noise1.wav");
v1 = resample(v1_20k, Fs, fs);
if size(v1, 2) > 1, v1 = mean(v1, 2); end
v1 = v1(:);

if length(v1) >= N_campioni
    all_noises(:,2) = v1(1:N_campioni);
else
    all_noises(:,2) = [v1; zeros(N_campioni - length(v1), 1)];
end

% --- 3. Medical Beep 2 ---
[v2_20k, fs] = audioread("rumori\Heart_monitor_beeping\medical_noise2.wav");
v2 = resample(v2_20k, Fs, fs);
if size(v2, 2) > 1, v2 = mean(v2, 2); end
v2 = v2(:);

if length(v2) >= N_campioni
    all_noises(:,3) = v2(1:N_campioni);
else
    all_noises(:,3) = [v2; zeros(N_campioni - length(v2), 1)];
end

% --- Parametri STFT ---
win_length = 512;
win = hamming(win_length); 
overlap = round(win_length * 0.75); 
nfft = 512;

% Calcolo STFT coerente su N_campioni = 60000 (Genera matrice K x 469)
V_home_alarm_temp = stft(all_noises(:,1), Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
V_beep1_temp      = stft(all_noises(:,2), Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);
V_beep2_temp      = stft(all_noises(:,3), Fs, "Window", win, "FFTLength", nfft, "OverlapLength", overlap);

% Limitiamo colonne come quelle degli speech
V_home_alarm = V_home_alarm_temp(:,1:449); 
V_beep1 = V_beep1_temp(:,1:449); 
V_beep2 = V_beep2_temp(:,1:449);

end