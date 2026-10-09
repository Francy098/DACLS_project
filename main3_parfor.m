%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcolo metriche per grafici Fig.7 pag.14

% Il database considerato è composto da 10 audio: 5 femminili e 5 maschili.
% Questi vengono corrotti con 3 tipologie di rumore (medical beep 1, medical beep 2, home alarm) affinché lo SNR sia 10 dB. In totale 
% quindi otteniamo 30 audio rumorosi (10 per ogni tipo di rumore). 
% Il calcolo viene effettuato con la lunghezza del filtro L=3 variando il
% ritardo tau = [1,10]. Per ogni tau calcoliamo le prestazioni dell'algoritmo NLMS, FMLNLMS, BMLNLMS e quello combinato.
% Per ognuno dei 30 audio, salviamo il calcolo delle varie metriche per poi rappresentarle opportunatamente, tramite  "Grafici.m". 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Avvia il pool parallelo se non è già attivo
clear all; clc;

if isempty(gcp('nocreate'))
    parpool; % Utilizza tutti i core fisici disponibili
end

%% Inizializzazione parametri

L = 3; % tappi filtro
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
%Per la PESQ measure, ogni file .wav vengono riscritti ad ogni ciclo di tau


%% Lettura Speech rumorosi e PULITI  

[y_all_noisy_speeches, y_all_clean_speeches] = my_speech_matrix(Fs);

%% Lettura rumori 
% Questa funzione si occupa di leggere tutti i rumori usati per rendere
% rumorosi i vari speech (home_alarm, medical beep 1 e 2), limitarne la
% lunghezza dei campioni a 60000 (come gli speech) e di restituire i loro
% valori nel dominio della freq. tramite la STFT. 
% Questi V_home_alarm, V_beep1, V_beep2 ci serviranno poi nel Noise
% Detector.

[V_home_alarm, V_beep1, V_beep2] = my_three_noises(Fs);

%% INIZIO FILTRAGGIO

parfor speech_index = 1:30
%for speech_index = 1:1         %TEST: per debug       

    fprintf('Noisy Speech %d ...\n', speech_index);

    %Salvataggio segnale rumoroso corrente 
    y_noisy = y_all_noisy_speeches(:,speech_index); 
    
    %Calcolo STFT di un noisy speech alla volta
    %STFT -> otteniamo una Y(k,m). cioè un matrice 2D con k = indice freq e m = indice frame
    Y1_temp = stft(y_noisy,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap); %ogni Fs c'è un valore sul vettore
    Y1 = Y1_temp(:,1:449); % LIMITE 

    [K, M] = size(Y1); % K = indice di frequenza, M = numero totale di frame
    
    % Preallocazione matrici
    E_f = zeros(K, M); Z_f = zeros(K, M); Hw_f = zeros(L, K, M);
    E_b = zeros(K, M); Z_b = zeros(K, M); Hw_b = zeros(L, K, M);
    E_c = zeros(K, M); Z_c = zeros(K, M); Hw_c = zeros(L, K, M);
    
    % Ottieniamo la lunghezza esatta dell'istft per questo specifico segnale
    y1_temp = istft(Y1,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap);
    len_istft = length(y1_temp);

    % Creazione degraded.wav corrente (PESQ)
    degraded_file = sprintf('degraded_%d.wav', speech_index);
    degraded_path = fullfile('sounds', degraded_file);
    audiowrite(degraded_path, y_noisy(1:len_istft), Fs);
    
    % Preallocazione dinamica (per evitare mismatch di dimensioni tra segnali diversi)
    e_f = zeros(len_istft, 10);
    e_b = zeros(len_istft, 10);
    e_c = zeros(len_istft, 10);
    
    % Preallocazione metriche
    vsd_f = zeros(10,1); zita_nrH_f = zeros(10,1); PESQ_f = zeros(10,1); STOI_f = zeros(10,1);
    vsd_b = zeros(10,1); zita_nrH_b = zeros(10,1); PESQ_b = zeros(10,1); STOI_b = zeros(10,1);
    vsd_c = zeros(10,1); zita_nrH_c = zeros(10,1); PESQ_c = zeros(10,1); STOI_c = zeros(10,1);
    % Preallocazione metriche per il non processato
    vsd_np = zeros(10,1); zita_nrH_np = zeros(10,1); PESQ_np = zeros(10,1); STOI_np = zeros(10,1);
    % Creazione filtro Hw_np --> per il rappresentazione metriche del non processato
    % Siccome E = Y1, basta implementare un Hw_np pieno di zeri
    Hw_np = zeros(L, K, M); 
    Lnp_bank = zeros(K, M);

    % Determinazione deterministica dell'indice del segnale pulito di riferimento
    clean_index = mod(speech_index - 1, 10) + 1;
    reference_full = y_all_clean_speeches(:, clean_index);
    reference_stoi = reference_full(1:len_istft); % Taglia alla lunghezza coerente con l'istft

    % --- Generazione del file reference pulito temporaneo nella cartella 'sounds' (PESQ) ---
    clean_file = sprintf('clean_%d.wav', speech_index);
    clean_path = fullfile('sounds', clean_file);
    audiowrite(clean_path, reference_stoi, Fs);

    %%%%%%%
    % STFT dello speech PULITO corrente, quindi NON RUMOROSO. Ci servirà in my_sperformance_measure
    X_clean_temp = stft(y_all_clean_speeches(:,clean_index),Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap); 
    X_clean = X_clean_temp(:,1:449); % LIMITE

    %fprintf('Start filtering!\n');

    for tau = 1:10 %Campioni di ritardo
    
        %fprintf('tau=%d ...\n', tau);

        %Determinazione del tipo di rumore per noise detector
        if(speech_index <= 10)
            V_noise = V_home_alarm;
        elseif((speech_index > 10)&&(speech_index <= 20))
            V_noise = V_beep1;
        else
            V_noise = V_beep2;
        end

        %Adaptive Filtering 
        [E_f,Z_f,Hw_f, E_b,Z_b,Hw_b, E_c,Z_c,Hw_c, Lf_bank,Lb_bank,Lc_bank] = Adaptive_Filtering(tau,L,Y1,V_noise);
    
        e_f(:,tau) = istft(E_f,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap);
        e_b(:,tau) = istft(E_b,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap);
        e_c(:,tau) = istft(E_c,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap);
    
        % File filtrati univoci per il calcolo parallelo del PESQ
        path_f = sprintf('sounds/filtered_f_%d_%d.wav', speech_index, tau);
        path_b = sprintf('sounds/filtered_b_%d_%d.wav', speech_index, tau);
        path_c = sprintf('sounds/filtered_c_%d_%d.wav', speech_index, tau);
        
        % Controllo che il max valore sia <= 1, altrimenti data clipping con audiowrite
        max_f = max(abs(e_f(:,tau)));
        max_b = max(abs(e_b(:,tau)));
        max_c = max(abs(e_c(:,tau)));
        if(max_f > 1)   e_f(:,tau) = e_f(:,tau)/max_f; end %Normalizzazione (se necessaria)
        if(max_b > 1)   e_b(:,tau) = e_b(:,tau)/max_b; end
        if(max_c > 1)   e_c(:,tau) = e_c(:,tau)/max_c; end

        % Ci servono i .wav per il calcolo del PESQ
        audiowrite(path_f, e_f(:,tau), Fs);
        audiowrite(path_b, e_b(:,tau), Fs);
        audiowrite(path_c, e_c(:,tau), Fs);  

        %% Performance measures

        [zita_nrH_f(tau,1),vsd_f(tau,1), ~] = my_performance_measures(tau,L,Hw_f,V_noise,X_clean,Lf_bank); %Forward
        [zita_nrH_b(tau,1),vsd_b(tau,1), ~] = my_performance_measures(tau,L,Hw_b,V_noise,X_clean,Lb_bank); %Backward
        [zita_nrH_c(tau,1),vsd_c(tau,1), ~] = my_performance_measures(tau,L,Hw_c,V_noise,X_clean,Lc_bank); %Combined Approch
        [zita_nrH_np(tau,1),vsd_np(tau,1), ~] = my_performance_measures(tau,L,Hw_np,V_noise,X_clean,Lnp_bank); %Not processed
        % Escludiamo gli NRR che non ci servono qui
        
        %[zita_nrH_f(tau,1),vsd_f(tau,1), ~] = Adaptive_Performance_Measures(tau,L,Hw_f,V_noise,X_clean,Lf_bank); %Forward
        %[zita_nrH_b(tau,1),vsd_b(tau,1), ~] = Adaptive_Performance_Measures(tau,L,Hw_b,V_noise,X_clean,Lb_bank); %Backward
        %[zita_nrH_c(tau,1),vsd_c(tau,1), ~] = Adaptive_Performance_Measures(tau,L,Hw_c,V_noise,X_clean,Lc_bank); %Combined Approch
        %[zita_nrH_np(tau,1),vsd_np(tau,1), ~] = Adaptive_Performance_Measures(tau,L,Hw_np,V_noise,X_clean,Lnp_bank); %Not processed
        
 
        % --- PESQ measure
        
        % Nomi dei file relativi (da passare a pesq2_mtlb)
        file_name_f = sprintf('filtered_f_%d_%d.wav', speech_index, tau);
        file_name_b = sprintf('filtered_b_%d_%d.wav', speech_index, tau);
        file_name_c = sprintf('filtered_c_%d_%d.wav', speech_index, tau);
        %audiowrite('sounds\degraded.wav',y_all_noisy_speeches(:,speech_index),Fs); %Creazione degraded.wav corrente
        PESQ_f(tau,1) = pesq2_mtlb( clean_file, file_name_f, Fs, 'wb', binary, pathaudio ); %full-band
        PESQ_b(tau,1) = pesq2_mtlb( clean_file, file_name_b, Fs, 'wb', binary, pathaudio ); %full-band
        PESQ_c(tau,1) = pesq2_mtlb( clean_file, file_name_c, Fs, 'wb', binary, pathaudio ); %full-band
        PESQ_np(tau)  = pesq2_mtlb( clean_file, degraded_file, Fs, 'wb', binary, pathaudio); % unprocessed
        
        %STOI measure
        STOI_f(tau) = stoi(e_f(:,tau), reference_stoi, Fs);
        STOI_b(tau) = stoi(e_b(:,tau), reference_stoi, Fs);
        STOI_c(tau) = stoi(e_c(:,tau), reference_stoi, Fs);
        STOI_np(tau) = stoi(y_noisy(1:len_istft), reference_stoi, Fs); % unprocessed
        
        % Pulizia opzionale sul disco per non accumulare giga di file temporanei
        delete(path_f); 
        delete(path_b);
        delete(path_c); 

    end

    %fprintf('End filtering!\n');
    
    delete(degraded_path); % Cancella il degraded temporaneo a fine ciclo
    delete(clean_path); % Eliminazione del file pulito temporaneo a fine ciclo

    % Salva solo le metriche di performance ignorando il filtro pesante Hw
    path_save = sprintf('results/Metriche_%d.mat',speech_index);
    % Nota: all'interno di parfor è consigliabile usare funzioni di salvataggio dedicate 
    % o salvare direttamente la struttura per evitare avvisi di performance.
    parsave(path_save, STOI_f, STOI_b, STOI_c, STOI_np, PESQ_f, PESQ_b, PESQ_c, PESQ_np, ...
                     zita_nrH_f, zita_nrH_b, zita_nrH_c, zita_nrH_np, ...
                     vsd_f, vsd_b, vsd_c, vsd_np, e_f, e_b, e_c);
   
    fprintf('Saved in Metriche_%d.mat\n',speech_index);

end

fprintf('All 30 noisy speeches successfully processed!\n');
