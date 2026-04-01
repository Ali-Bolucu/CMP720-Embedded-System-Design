classdef Visualizer
    methods (Static)
        
        %% 1. Pareto Cephesi Grafiği (MCMADisper uyarlaması)
        function plotPareto(AllObjectives, outputDir)
            % Arka planda figür oluştur (Ekranda pop-up çıkmasını engeller)
            fig = figure('Name', 'Pareto Front', 'Visible', 'off');
            
            % Obj(:, 1) -> Latency, Obj(:, 2) -> Fault Tolerance
            scatter(AllObjectives(:, 1), AllObjectives(:, 2), 50, 'k', 'filled');
            
            xlabel('Latency', 'FontSize', 14);
            ylabel('Fault Tolerance', 'FontSize', 14);
            title('Pareto Front (NSGA-II)', 'FontSize', 14);
            grid on;
            
            % Çıktıyı kaydet
            saveas(fig, fullfile(outputDir, 'Pareto_Front.png'));
            close(fig);
        end
        
        %% 2. NoC Yerleşim Grafiği (MOstraNoC uyarlaması)
        function plotNoC(bestDecision, MegaGraph, outputDir)
            fig = figure('Name', 'NoC Mapping', 'Visible', 'off');
            
            % Parametreleri MegaGraph'tan çek
            num_linhas = cell2mat(MegaGraph(2));
            num_colunas = cell2mat(MegaGraph(3));
            WTask = cell2mat(MegaGraph(8));
            num_App = length(WTask);
            
            % Uygulamaların görev sınırlarını belirle
            tarefa_app = zeros(2, num_App);
            for j = 1:num_App
                if j == 1
                    tarefa_app(1,j) = 1;
                    tarefa_app(2,j) = WTask(j);
                else
                    tarefa_app(1,j) = WTask(j-1) + 1;
                    tarefa_app(2,j) = WTask(j);
                end
            end
            
            % Görsel ayarları
            Cores = [0 .5 .5; 0 .35 .850; .33 .65 .65; .45 .55 .65; 0.5 0 0.5];
            width = 2;       
            height = 2;      
            
            % Görevleri Çekirdeklere (Cores) Dağıt
            PosiNoc = zeros((num_linhas * num_colunas), 1);
            for s = 1:(num_linhas * num_colunas)
                 if ismember(s, bestDecision)
                     Var1 = find(bestDecision(1,:) == s);
                     PosiNoc(s, 1:length(Var1)) = Var1;
                 end
            end
            
            % Matrisi NoC boyutlarına göre şekillendir
            max_tasks_per_core = size(PosiNoc, 2);
            numeros = reshape(PosiNoc, num_linhas, num_colunas, max_tasks_per_core);
            numeros = permute(numeros, [2, 1, 3]);
            
            % Izgarayı çiz
            for i = 1:num_linhas
                for j = 1:num_colunas
                    x = width*(j - 1);
                    y = height*(num_linhas - i);
                    
                    Corquadrado = [1 1 1]; % Boş çekirdek rengi (Beyaz)
                    
                    % Çekirdekteki ilk göreve göre renk belirle
                    if numeros(i, j, 1) ~= 0
                        first_task = numeros(i, j, 1);
                        for l = 1 : num_App
                            if first_task >= tarefa_app(1, l) && first_task <= tarefa_app(2, l)
                                % Renk paletinde sınır aşımını önlemek için mod al
                                color_idx = mod(l-1, size(Cores,1)) + 1; 
                                Corquadrado = Cores(color_idx, :);
                            end
                        end
                    end
                    
                    rectangle('Position', [x y width height], 'FaceColor', Corquadrado, 'LineWidth', 2);
                    
                    % Çekirdeğe atanan tüm görevleri yazdır
                    tasks_on_core = nonzeros(numeros(i, j, :));
                    if ~isempty(tasks_on_core)
                        task_str = sprintf('%d ', tasks_on_core');
                        text(x + width/2, y + height/2, task_str, 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
                    end
                end
            end
            
            axis([-1 num_colunas*width+1 -1 num_linhas*height+1]);
            axis off;
            
            saveas(fig, fullfile(outputDir, 'NoC_Mapping.png'));
            close(fig);
        end
        
    end
end