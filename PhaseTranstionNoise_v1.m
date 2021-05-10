% Gera o Diagrama de Transição de Fase
% Utilizado no artigo submetido ao SBRT 2020

% Fez uso do código:

% l1eq_example.m
%
% Test out l1eq code (l1 minimization with equality constraints).
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

% put key subdirectories in path if not already there
path(path, './Optimization');
path(path, './Data');

% Tamanho do sinal
n = 80;

% Quantidade de iterações

N_iteracoes = 1000;

% Matriz com as taxas de sucesso
sucesso = zeros(n, n);

% SNR (em dB)
SNR =80;


% Caso queira um diagram baseado na Matriz Gaussiana descomente e linha abaixo
%A_TMP = randn(n,n);

      
% Caso queira um diagrama com uma matriz DCT descomente a linha abaixo
A_TMP = dctmtx(n);

for m = 1:n % Varia a quantidade de medidas de 1 até n   
    for k = 1:1:m % Varia a esparsidade do sinal de 1 até a quantide de medidas
            erro_acum = 0;
            for it = 1:N_iteracoes % Quantidede de iterações
                % Gera o sinal com valores +/- 1
                x = zeros(n,1);
                T = randi([1,k]);
                q = randperm(T);
                x(q(1:T)) = sign(randn(T,1));
                
                sigma_2 = norm(x)*10^(-SNR/10); % Variância do ruído que afeta o sinal
                sigma0_2 = norm(x)*10^(-SNR/10); % Variância do ruído que afeta a medição
                % Adicionando ruido ao sinal
                x = x + sqrt(sigma_2).*randn(n, 1);
              
                rows = randperm(n);
                A = A_TMP(rows(1:m),:);

                % Medidas com a possibilidade de adicionar ruído se
                % sigma0_2 for diferente de zero.
                y = A*x + sqrt(sigma0_2).*randn(m, 1);
                
                % Chute inicial = min energy
                x0 =     A'*y;
                
                % Resolve a minização
                tic
                xp = l1eq_pd(x0, A, [], y, 1e-3);
                toc
                                
                erro = norm(x-xp)/norm(x); % Calcula o erro
                
                % large scale
                % Afun = @(z) A*z;
                % Atfun = @(z) A'*z;
                % tic
                % xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);
                % toc
                
                if (erro < 1e-3) % Conta os sucessos
                    sucesso(n - k + 1, m) = sucesso(n - k + 1, m) + 1;
                end
            end % Fim das iteracoes para um dado k
    end %fim variacao da esparsidade (k).
end % Fim variacao do numero de medidas (m)

% Salva a matriz com os resultados
save('n80SNR80Sinalheatmap.mat','n','sucesso','N_iteracoes')

%% Desenho do heatmap como uma região

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

 %% Desenho das curvas no heatmap
 
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