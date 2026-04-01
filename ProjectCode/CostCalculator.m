classdef CostCalculator
    methods(Static)
        
        %% Fault Tolerance
        function ftCost = FaultToleranceCost(PopDec, numRows, numCols, DistTab)
            
            nSolutions = size(PopDec, 1);
            nCores = numRows * numCols;
            ftCost = zeros(nSolutions, 1);
            allCores = 1:nCores;
            
            for i = 1:nSolutions
                % Dolu ve boş çekirdekleri bul
                busyCores = unique(PopDec(i, :));
                
                % Hata Ayıklama: busyCores içinde 0 veya negatif değer var mı?
                if any(busyCores <= 0)
                    warning('FaultToleranceCost: Solution %d contains invalid core assignments (<= 0).', i);
                    busyCores = busyCores(busyCores > 0); % Sadece geçerli olanları tut
                end
                
                idleCores = setdiff(allCores, busyCores);
                
                % Eğer boş (yedek) çekirdek kalmamışsa algoritmaya çok büyük bir ceza ver
                if isempty(idleCores)
                    ftCost(i) = 99999; % Maksimum ceza puanı
                    continue;
                end
                
                % --- VEKTÖREL SİHİR BURADA ---
                % DistTab tablosundan, sadece DOLU çekirdeklerin BOŞ çekirdeklere 
                % olan tüm mesafelerini tek hamlede çek.
                % Sonra min(..., [], 2) ile her dolu çekirdek için EN KISA mesafeyi bul.
                minDistances = min(DistTab(busyCores, idleCores), [], 2);
                
                % Ftol = distance - 1 (1 hücre uzaklık = 0 maliyet)
                penalties = minDistances - 1;
                
                % Toplam FT maliyetini kaydet
                ftCost(i) = sum(penalties);
            end
        end

        %% Latency Cost
        function latencyCost = LatencyCost(Population, DistanceMatrix, SourceTasks, TargetTasks, EdgeWeights)
            
            nSolutions = size(Population, 1);
            nTasks = size(Population, 2);
            nCores = size(DistanceMatrix, 1);
            latencyCost = zeros(nSolutions, 1);
            
            % Validate input task indices
            if any(SourceTasks < 1 | SourceTasks > nTasks) || any(TargetTasks < 1 | TargetTasks > nTasks)
                % Invalid task indices - assign high penalty
                latencyCost(:) = 999999;
                warning('Invalid task indices in LatencyCost: Source [%d,%d], Target [%d,%d], nTasks=%d', ...
                    min(SourceTasks), max(SourceTasks), min(TargetTasks), max(TargetTasks), nTasks);
                return;
            end
            
            % Validate edge weights length matches task pairs
            if length(EdgeWeights) ~= length(SourceTasks)
                latencyCost(:) = 999999;
                warning('EdgeWeights length (%d) does not match edge count (%d)', length(EdgeWeights), length(SourceTasks));
                return;
            end
            
            for i = 1:nSolutions
                % Bu çözümdeki TÜM kaynak ve hedef görevlerin atandığı çekirdekleri TEK HAMLEDE bul
                sourceCores = Population(i, SourceTasks);
                targetCores = Population(i, TargetTasks);
                
                % Geçerliliği kontrol et: Tüm çekirdek indisleri [1, nCores] aralığında olmalı
                if any(sourceCores < 1 | sourceCores > nCores) || any(targetCores < 1 | targetCores > nCores)
                    % Geçersiz çekirdek ataması - çok büyük ceza ver
                    latencyCost(i) = 999999;
                    continue;
                end
                
                % Bu çekirdekler arasındaki tüm mesafeleri Dist_Tab içinden TEK HAMLEDE çek
                linearIndices = sub2ind(size(DistanceMatrix), sourceCores, targetCores);
                distances = DistanceMatrix(linearIndices);
                
                % Latency = Edge weight + Distance contribution (Sizin formülünüz)
                % Not: Literatür standardı için burayı (EdgeWeights .* distances) yapabilirsiniz.
                edgeLatencies = EdgeWeights + (distances .* 10);
                
                % Tüm kenarların maliyetini topla ve normalize et
                latencyCost(i) = sum(edgeLatencies) / 1000; 
            end
            
        end
        
    end
end