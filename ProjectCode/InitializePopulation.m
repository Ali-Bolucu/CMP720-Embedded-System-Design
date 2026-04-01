%% InitializePopulation.m


function PopOut = InitializePopulation(Npop, TotalTasks, CoresRow, CoresColumn, AppLengths)
            
            % NoC boyutlarını hesapla
            numRows = CoresRow;
            numCols = CoresColumn;
            numCores = CoresRow * CoresColumn;
            
            % Bu şekilde başlangıç ve bitiş noktaları tutuldu
            % endApp = [16, 28, 53, 73]
            % startApp = [1, 17, 29, 54]
            endApp = AppLengths;                     
            startApp = [1, endApp(1:end-1) + 1];
            
            % Manhattan mesafe matrisini oluştur
            [LN, CL] = ind2sub([numRows, numCols], 1:numCores);  % Her node'u bir X,Y'ye atadı
            Pos_Tab = [LN' CL'];                                 % Matrixe dönüştürdü, yukarda 2 vektör vardı
            Dist_Tab = pdist2(Pos_Tab, Pos_Tab, 'cityblock');    % Her 2 node arasındaki mesafe
            
            PopOut = zeros(Npop, TotalTasks); % Çözüm kümesi: Npop x 73 (100 çözüm, her biri 73 görev için çekirdek ataması)
            
            for i = 1 : Npop
                gridNoC = zeros(numRows, numCols);
                randAppOrder = randperm(length(AppLengths)); % Uygulamaları hangi sırayla yerleştireceğiz
                
                for k = 1 : length(AppLengths)
                    appIndex = randAppOrder(k); % Uygulamaların yerleştirme sırasını karıştırmak için
                    
                    % -------------------------------------------------------------
                    % BURASI SPEKTRAL KÜMELEMENİN GELECEĞİ YERDİR!
                    % Şimdilik eski kodun yaptığı gibi görevleri rastgele diziyoruz:
                    taskRange = startApp(appIndex) : endApp(appIndex);
                    shuffledTasks = taskRange(randperm(length(taskRange)));
                    % İleride burayı: "shuffledTasks = helper.SpectralClustering(...)" 
                    % şeklinde değiştireceğiz.
                    % -------------------------------------------------------------
                    
                    for m = 1 : length(shuffledTasks)
                        currentTask = shuffledTasks(m);
                        
                        if m == 1
                            % İlk görev için rastgele bir (X,Y) bul
                            randX = randi(numRows, 1, 1);
                            randY = randi(numCols, 1, 1);
                            
                            % Eğer doluysa, en yakın boş yeri bul
                            if gridNoC(randX, randY) ~= 0
                                [randX, randY] = findNearestFreeCore(gridNoC, randX, randY, Dist_Tab, numRows, numCols);
                            end
                            gridNoC(randX, randY) = currentTask;
                        else
                            % Diğer görevleri bir önceki görevin yakınına yerleştir
                            [randX, randY] = findNearestFreeCore(gridNoC, randX, randY, Dist_Tab, numRows, numCols);
                            gridNoC(randX, randY) = currentTask;
                        end
                    end
                end
                
                % 2D NoC matrisini 1D Kromozoma (Diziye) çevir
                % PopOut matrisinin her satırı [Görev 1'in Çekirdeği, Görev 2'nin Çekirdeği, ...] şeklinde olmalı
                flatNoC = reshape(gridNoC', 1, numRows * numCols);
                chromosome = zeros(1, TotalTasks);
                for coreID = 1:(numRows * numCols)
                    taskID = flatNoC(coreID);
                    if taskID > 0
                        chromosome(taskID) = coreID;
                    end
                end
                PopOut(i,:) = chromosome;
            end
        end
        
    %% 2. En Yakın Boş Çekirdeği Bulma Fonksiyonu
        function [nextX, nextY] = findNearestFreeCore(gridNoC, currX, currY, DistTab, numRows, numCols)
            % DÜZELTME 2: gridNoC yanındaki transpoze (') işaretini sildik
            flatNoC = reshape(gridNoC, 1, numRows * numCols);
            currentCoreID = sub2ind([numRows, numCols], currX, currY); 
            
            maxDist = max(DistTab(currentCoreID, :));
            
            for dist = 1 : maxDist
                % Belli bir uzaklıktaki tüm çekirdekleri bul
                coresAtDist = find(DistTab(currentCoreID, :) == dist);
                % Bu çekirdeklerden boş olanları (0 olanları) kesiştir
                freeCores = intersect(coresAtDist, find(flatNoC == 0));
                
                if ~isempty(freeCores)
                    % Eğer birden fazla boş yer varsa rastgele birini seç
                    chosenCoreID = freeCores(randi(length(freeCores), 1, 1));
                    
                    [nextX, nextY] = ind2sub([numRows, numCols], chosenCoreID);
                    return;
                end
            end
            
            % Eğer NoC tamamen dolduysa
            nextX = 1; nextY = 1; 
        end
        