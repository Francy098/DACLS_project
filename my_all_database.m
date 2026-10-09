%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Questa funzione è pensata per fornire in uscita due matrici:
% - all_speeches --> dove vi sono tutti i 20 speech puliti
% - all_noises ----> dove ci sono tutti i 26 rumori considerati
% Tutti loro sono campionati alla freq Fs = 16kHz e cosideriamo per tutti
% 60000 campioni.

function [my_all_speeches, my_all_noises] = my_all_database(Fs)

N_campioni = 60000;

% Inizializzazione matrici
my_all_speeches = zeros(N_campioni, 20);
my_all_noises = zeros(N_campioni, 26);


%% Inserimento voci maschili in my_all_speeches
%
for i = 1:10
    path_audio = sprintf("our_TIMIT_database/male_%d.wav",i);
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
    my_all_speeches(:,i) = y_temp2;
end   

%% Inserimento voci femminili in my_all_speeches
%
for i = 1:10
    path_audio = sprintf("our_TIMIT_database/female_%d.wav",i);
    [y,fs] = audioread(path_audio);
    
    % Ricampionamento del segnale da a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else %altrimenti Zero padding
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; 
    end

    % La salviamo in y home alarm
    my_all_speeches(:,i + 10) = y_temp2;
end


%% Inserimento RUMORI in my_all_noises
%
for i = 1:26
    path_audio = sprintf("rumori/noise_%d.wav",i);
    [y,fs] = audioread(path_audio);
    
    % Ricampionamento del segnale da a 16 kHz (Fs)
    y_temp = resample(y, Fs, fs);
    y_temp = y_temp(:); % Assicura che sia un vettore colonna [N x 1]
    
    % Uniformiamo la lunghezza a N_campioni
    if length(y_temp) >= N_campioni
        y_temp2 = y_temp(1:N_campioni); %troppo lungo? Allora tagliamo
    else %altrimenti Zero padding
        y_temp2 = [y_temp; zeros(N_campioni - length(y_temp), 1)]; 
    end

    % La salviamo in y home alarm
    my_all_noises(:,i) = y_temp2;
end

end