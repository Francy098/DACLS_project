%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcolo metriche per grafici Fig.8 pag.15
%
% Come databese utilizzeremo solo gli speech corrotti con il rumore
% hospital beeping 1 (Beep1). Verranno processati e si considereranno solo
% le metriche e prestazioni dell'approccio combinato per L = {3,5,8,10} con
% tau = [1:10] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inizializzazione parametri
clear all; clc;

L_vec = [3;5;8;10];
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

% Prendiamo solo gli speech corrotti dal Beep 1
y_our_noisy_speeches = y_all_noisy_speeches(:, 11:20);

%% Lettura rumore interessato Beep1 
[ ~, V_noise, ~] = my_three_noises(Fs);

%% Preallocazione metriche 
vsd_c = zeros(10,4); % 4 perché gli L considerati sono 4
zita_nrH_c = zeros(10,4); 
PESQ_c = zeros(10,4); 
STOI_c = zeros(10,4);
PESQ_np = zeros(10,4); %unprocessed
STOI_np = zeros(10,4);


%% INIZIO FILTRAGGIO

fprintf('Start filtering!\n');

for idx_L = 1:4 % 4 perché L_vec è lungo 4

    L_curr= L_vec(idx_L);
    fprintf('Elaborazione per L = %d...\n', L_curr);
    
    % Preallocazione metriche da usare per ogni L
    vsd_c_Lcurr = zeros(10,10); 
    zita_nrH_c_Lcurr = zeros(10,10); 
    PESQ_c_Lcurr = zeros(10,10); 
    STOI_c_Lcurr = zeros(10,10);
    PESQ_np_Lcurr = zeros(10,10); %unprocessed
    STOI_np_Lcurr = zeros(10,10);

    for speech_index = 1:10
           
        fprintf('Noisy Speech %d ...\n', speech_index);
    
        %Salvataggio segnale rumoroso corrente 
        y_noisy = y_our_noisy_speeches(:,speech_index); 
        
        %Calcolo STFT di un noisy speech alla volta
        %STFT -> otteniamo una Y(k,m). cioè un matrice 2D con k = indice freq e m = indice frame
        Y1_temp = stft(y_noisy,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap); %ogni Fs c'è un valore sul vettore
        Y1 = Y1_temp(:,1:449); % LIMITAZIONE per far combaciare dimensioni 
        [K, M] = size(Y1); % K = indice di frequenza, M = numero totale di frame
        
        % Preallocazione matrici
        E_c = zeros(K, M); Z_c = zeros(K, M); Hw_c = zeros(L_curr, K, M);
    
        % Creazione degraded.wav corrente per PESQ measure
        degraded_file = sprintf('degraded_%d.wav', speech_index);
        degraded_path = fullfile('sounds', degraded_file);
        audiowrite(degraded_path, y_noisy, Fs);
        
        % Ottieniamo la lunghezza esatta dell'istft per questo specifico segnale
        y1_temp = istft(Y1,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap);
        len_istft = length(y1_temp);

        % Preallocazione dinamica (per evitare mismatch di dimensioni tra segnali diversi)
        e_c = zeros(len_istft, 10);
        
        % Determinazione deterministica dell'indice del segnale pulito di riferimento
        reference_full = y_all_clean_speeches(:, speech_index);
        reference_stoi = reference_full(1:len_istft); % Taglia alla lunghezza coerente con l'istft
    
        % --- Generazione del file reference pulito temporaneo nella cartella 'sounds' (PESQ) ---
        clean_file = sprintf('clean_%d.wav', speech_index);
        clean_path = fullfile('sounds', clean_file);
        audiowrite(clean_path, reference_stoi, Fs);

        %%%%%%%
        % STFT dello speech PULITO corrente, quindi NON RUMOROSO. Ci servirà in my_sperformance_measure
        X_clean = stft(y_all_clean_speeches(:,speech_index),Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap); 
    
    
        %fprintf('Start filtering!\n');
    
        for tau = 1:10 %Campioni di ritardo
        
            %fprintf('tau=%d ...\n', tau);
    
            %Adaptive Filtering 
            [~,~,~, ~,~,~, E_c,Z_c,Hw_c, Lf_bank,Lb_bank,Lc_bank] = Adaptive_Filtering(tau, L_curr, Y1, V_noise);
        
            e_c(:,tau) = istft(E_c,Fs,"Window",win,"FFTLength",nfft,"OverlapLength",overlap);
        
            % File filtrati univoci per il calcolo parallelo del PESQ
            path_c = sprintf('sounds/filtered_c_%d_%d.wav', speech_index, tau);
            
    
            % Controllo che il max valore sia <= 1, altrimenti data clipping con audiowrite
            max_c = max(abs(e_c(:,tau)));
            if(max_c > 1)   e_c(:,tau) = e_c(:,tau)/max_c; end
    
            % Ci servono i .wav per il calcolo del PESQ
            audiowrite(path_c, e_c(:,tau), Fs);  
    
            %% Performance measures

           [zita_nrH_c_Lcurr(tau,speech_index), vsd_c_Lcurr(tau,speech_index), ~] = my_performance_measures(tau,L_curr,Hw_c,V_noise,X_clean,Lc_bank); %Combined Approch
            % Escludiamo i gamma_x e gamma_v e NRR che non ci servono qui
    
            % PESQ measure
            % Nomi dei file relativi (da passare a pesq2_mtlb)
            file_name_c = sprintf('filtered_c_%d_%d.wav', speech_index, tau);
            PESQ_c_Lcurr(tau,speech_index)  = pesq2_mtlb(clean_file, file_name_c, Fs, 'wb', binary, pathaudio ); %full-band
            PESQ_np_Lcurr(tau,speech_index) = pesq2_mtlb(clean_file, degraded_file, Fs, 'wb', binary, pathaudio ); %full-band

            %STOI measure
            STOI_c_Lcurr(tau,speech_index) = stoi(e_c(:,tau), reference_stoi, Fs);
            STOI_np_Lcurr(tau,speech_index) = stoi(y_noisy(1:len_istft), reference_stoi, Fs);
            
            % Pulizia opzionale sul disco per non accumulare giga di file temporanei
            delete(path_c); 
            
        end

        %fprintf('End filtering!\n');
        delete(degraded_path); % Cancella il degraded temporaneo a fine ciclo del segnale
        delete(clean_path); 
     end
        
     % Calcolo della corrispetiva metrica _Lcurr
     vsd_c(:,idx_L) = mean(vsd_c_Lcurr,2);
     zita_nrH_c(:,idx_L) = mean(zita_nrH_c_Lcurr,2); 
     PESQ_c(:,idx_L) = mean(PESQ_c_Lcurr,2); 
     PESQ_np(:,idx_L) = mean(PESQ_np_Lcurr,2);
     STOI_c(:,idx_L) = mean(STOI_c_Lcurr,2);
     STOI_np(:,idx_L) = mean(STOI_np_Lcurr,2);
end

fprintf('All noisy speeches successfully processed for all L values!\n');

%% Salva solo le metriche di performance ignorando il filtro pesante Hw
   save('results/Metriche_Adaptive_varius_L.mat', ...
        'vsd_c', 'zita_nrH_c', 'PESQ_c','PESQ_np','STOI_np', 'STOI_c');

fprintf('Saved in Metriche_Adaptive_VariusL.mat\n');

%% Grafici

tau_x_axis = 1:10;

load('results/Metriche_Adaptive_varius_L.mat');

% ----------- Inizializzazione metriche per il non processato -------------
% Come visto in precedenza, nel main3_parfor.m e con Grafici.m, il non
% processato presenta un vsd = 0 e con uno zita_dB = 0 dB. Quindi creiamo
% due vettori affinché vengano replicate queste condizioni

%vsd_np = zeros(1, 10); ---- 10log10(vsd_np) NON DEFINITO
zita_np_dB = zeros(1, 10);

%Normalizzazione vsd
vsd_c_norm = vsd_c(:,:)/max(vsd_c, [] ,'all');

% vsd_c in dB
vsd_c_dB = 10*log10(vsd_c_norm);

% % Fig. 8(a) Distortion Index vsd
figure;
subplot(2,2,1); %posizionamento dei grafici
plot(tau_x_axis,vsd_c_dB(:,1),'r-o',tau_x_axis,vsd_c_dB(:,2),'y-s',tau_x_axis,vsd_c_dB(:,3),'b-d',tau_x_axis,vsd_c_dB(:,4),'v-','LineWidth', 1.2);
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('L= 3','L= 5','L= 8','L= 10',Location= 'best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("v_{sd}[dB]",'FontWeight','bold');
title('Fig. 8(a) Distortion Index (v_{sd})');

% Portiamo tutto in dB
zita_nrH_c_dB = 10*log10(abs(zita_nrH_c));

% Fig. 8(b) Noise Reduction Factor zita_nr (in dB)
subplot(2,2,2); %posizionamento dei grafici
plot(tau_x_axis,zita_nrH_c_dB(:,1),'r-o',tau_x_axis,zita_nrH_c_dB(:,2),'y-s',tau_x_axis,zita_nrH_c_dB(:,3),'b-d',tau_x_axis,zita_nrH_c_dB(:,4),'v-',tau_x_axis,zita_np_dB,'g','LineWidth', 1.2);
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('L= 3','L= 5','L= 8','L= 10','unprocessed',Location= 'best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("\xi_{nr}[dB]",'FontWeight','bold');
title('Fig. 8(b) Noise Reduction Factor (\xi_{nr})');


%% Fig. 8(c) PESQ
%figure;
subplot(2,2,3); %posizionamento dei grafici
plot(tau_x_axis,PESQ_c(:,1),'r-o',tau_x_axis,PESQ_c(:,2),'y-s',tau_x_axis,PESQ_c(:,3),'b-d',tau_x_axis,PESQ_c(:,4),'v-',tau_x_axis,PESQ_np,'g','LineWidth', 1.2);
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('L= 3','L= 5','L= 8','L= 10','unprocessed',Location= 'best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("PESQ",'FontWeight','bold');
title('Fig. 8(c) PESQ');

% Fig. 8(d) STOI
%figure;
subplot(2,2,4); %posizionamento dei grafici
plot(tau_x_axis,STOI_c(:,1),'r-o',tau_x_axis,STOI_c(:,2),'y-s',tau_x_axis,STOI_c(:,3),'b-d',tau_x_axis,STOI_c(:,4),'v-',tau_x_axis,STOI_np,'g','LineWidth', 1.2);
grid on;
%xlim([0 10]); %per limitare asse x
%ylim([0 0.014]); %per limitare asse y
legend('L= 3','L= 5','L= 8','L= 10','unprocessed',Location= 'best'); %mette la legenda in modo tale che non copre il grafico
xlabel("\tau",'FontWeight','bold'); %in grassetto
ylabel("STOI",'FontWeight','bold');
title('Fig. 8(d) STOI');