function [y_all_noisy_speeches, y_all_clean_speeches] = my_speech_matrix(Fs)

N_campioni = 60000; % Massimo numero dei campioni che consideriamo
y_home_alarm = zeros(N_campioni, 10);
y_beep1 = zeros(N_campioni, 10);
y_beep2 = zeros(N_campioni, 10);

%% Creazione y_home_alarm

% Campioniamo le voci maschili
for i = 1:5
    path_audio = sprintf("sounds/home_alarm/noisy_male_%d_SNR10dB.wav",i);
    [y,fs] = audioread(path_audio);
    
    % Ricampionamento del segnale da a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y home alarm
    y_home_alarm(:,i) = y_temp2;
    
end

% Campioniamo le voci femminili
for i = 1:5
    path_audio = sprintf("sounds/home_alarm/noisy_female_%d_SNR10dB.wav",i);
    [y,fs] = audioread(path_audio);
                                   
    % Ricampionamento del segnale a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]

    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y_home_alarm, partendo da 6 perché da 1 a 5 ci sono le voci maschili
    y_home_alarm(:,i+5) = y_temp2;
end

%% Creazione y_beep1

% Campioniamo le voci maschili
for i = 1:5
    path_audio = sprintf("sounds/medical_beep1/noisy_male_%d_SNR10dB.wav",i);
    [y,fs] = audioread(path_audio);
                                    
    % Ricampionamento del segnale a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y home alarm
    y_beep1(:,i) = y_temp2;
end

% Campioniamo le voci femminili
for i = 1:5
    path_audio = sprintf("sounds/medical_beep1/noisy_female_%d_SNR10dB.wav",i);
    [y,fs] = audioread(path_audio); %Prendiamo max numero di campioni = 32313
                                   % fs=20kHz di default
    
    % Ricampionamento del segnale da 20 kHz (fs) a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]

    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y_home_alarm, partendo da 6 perché da 1 a 5 ci sono le voci maschili
    y_beep1(:,i+5) = y_temp2;
end

%% Creazione y_beep2

% Campioniamo le voci maschili
for i = 1:5
    path_audio = sprintf("sounds/medical_beep2/noisy_male_%d_SNR10dB.wav",i);
    [y,fs] = audioread(path_audio); %Prendiamo max numero di campioni = 32313
                                   % fs=20kHz di default
    
    % Ricampionamento del segnale da 20 kHz (fs) a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]

    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y home alarm
    y_beep2(:,i) = y_temp2;
end

% Campioniamo le voci femminili
for i = 1:5
    path_audio = sprintf("sounds/medical_beep2/noisy_female_%d_SNR10dB.wav",i);
    [y,fs] = audioread(path_audio); %Prendiamo max numero di campioni = 32313
                                   % fs=20kHz di default
    
    % Ricampionamento del segnale da 20 kHz (fs) a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y_home_alarm, partendo da 6 perché da 1 a 5 ci sono le voci maschili
    y_beep2(:,i+5) = y_temp2;
end

% Creazione di un'unica matrice in cui inseriamo tutti i noisy speech
y_all_noisy_speeches = [y_home_alarm, y_beep1, y_beep2];
%se al posto della ,(virgola) mettessimo ; verrebbero messi in riga


%% Creazione y_all_clean_speeches

%Per salvataggio degli speech puliti
y_all_clean_speeches = zeros(N_campioni, 10); %perché abbiamo 10 voci diverse

% Campioniamo le voci maschili
for i = 1:5
    path_audio = sprintf("our_TIMIT_database/male_%d.wav",i);
    [y,fs] = audioread(path_audio); %Prendiamo max numero di campioni = 32313
                                   % fs=20kHz di default
    
    % Ricampionamento del segnale da 20 kHz (fs) a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y home alarm
    y_all_clean_speeches(:,i) = y_temp2;
end

% Campioniamo le voci femminili
for i = 1:5
    path_audio = sprintf("our_TIMIT_database/female_%d.wav",i);
    [y,fs] = audioread(path_audio); %Prendiamo max numero di campioni = 32313
                                   % fs=20kHz di default
    
    % Ricampionamento del segnale da 20 kHz (fs) a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; %altrimenti Zero padding
    end

    % La salviamo in y_home_alarm, partendo da 6 perché da 1 a 5 ci sono le voci maschili
    y_all_clean_speeches(:,i+5) = y_temp2;
end


end