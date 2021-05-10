% Gera os Diagramas de Transição de Fase do artigo submetido ao SBRT 2020
% Utilização:
%  Cada sesção desenha um diagrama e está nomeado com o número da Figura no
%  texto. Para gerar uma figura basta descomentar a seção. 


 %% Figura 2

 % Tamanho do sinal usado no artigo
n = 80;

 load('n80heatmapm0.mat')
 
 h = heatmap(sucesso./N_iteracoes,'Xlabel','m/n','Ylabel','k/m');

% Ajustes nos rótulos
for x = 1:n
    if (mod(x, 20) ~= 0)
        CustomXLabels(x) = " ";
    else
        CustomXLabels(x) = string(x/n);
    end
end

for y = 1:n
    if (mod(y, 20) ~= 0)
        CustomYLabels(y) = " ";
    else
        CustomYLabels(y) = string((n-y)/n);
    end
end

 % Set the 'XDisplayLabels' property of the heatmap 
 % object 'h' to the custom x-axis tick labels
 h.XDisplayLabels = CustomXLabels;
 h.YDisplayLabels = CustomYLabels;
 
 grid off
 
  
 %% Figura 3
 
 % Tamanho do sinal usado no artigo
n = 80;

load('n80heatmapm0.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.4)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
  ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.9)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 %% Figura 4
 
 % Tamanho do sinal usado no artigo
n = 80;

load('n80heatmapm0.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 load('n80SNR80Sinalheatmap.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 %% Figura 5
 
 % Tamanho do sinal usado no artigo
n = 80;

load('n80heatmapm0.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 load('n80SNR80Medidaheatmap.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 %% Figura 6
 
 % Tamanho do sinal usado no artigo
n = 80;

load('n80heatmapm0.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 load('n80SNR80Medidaheatmap.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 
 load('n80SNR80Sinalheatmap.mat')
 
 ind = 1;
 for i=1:length(sucesso)
     j = 1;
     while(sucesso(i ,j)/N_iteracoes < 0.8)
         j = j + 1;
     end
     eixo_x(ind) = (j-1)/n;
     eixo_y(ind) = i/n;
     linha_1(ind) = sucesso(i, j-1)/N_iteracoes;
     ind = ind + 1;
 end

 for iter = 1:length(linha_1)
    grid on;
    hold on;          
    plot3(flip(eixo_x),eixo_y, linha_1);
    xlabel('m/n')
    ylabel('k/m')
    drawnow   
 end
 