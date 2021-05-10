% l1eq_example.m
%
% Test out l1eq code (l1 minimization with equality constraints).
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

clear all

% put key subdirectories in path if not already there
path(path, './Optimization');
path(path, './Data');

% To reproduce the example in the documentation, uncomment the 
% two lines below
%load RandomStates
%rand('state', rand_state);
%randn('state', randn_state);


% number of spikes in the signal
%T = randi([1,5]);
% number of observations to make
%m = 40;

% QUantidade de repetições
Iterations = 10;

% Vairação da quantidad de medidas
HorizontalDistance = 20;

% Limite para qual a probabilidade é aceitável. A partir desse 
% valor conidera-se que deve ser usado métodos gulosos.
Bound_S1 = 0.9;

% Limite entre o vetor recuperado e original considerado como 
% aceitável. Menor que esse valor, é considerado um acerto.
Bound_S2 = 1e-2;

% Inicialização da esparsidade.
k = 1;

% Tamanho do vetor.
n = 100;

% Inicialização da quantidade de medidas
m = HorizontalDistance;

% Taxa máxima de acerto. Pontos no qual a curva será desenhada.
SucessMap = zeros(n,n);

% Variável para guardar o valor da taxa de acerto.
acerto_anterior = 0;

% Relação Sinal-Ruído em dB
SNR = 75;

% Laço para determinação dos valores que comporão a curva.
while (m < n)
    acerto = 0;
    for it = 1:Iterations
        %Resolve o problema de minização
        % measurement matrix
        disp('Creating measurment matrix...');
        A = randn(m,n);
        A = orth(A')';
        disp('Done.');
        
        % random +/- 1 signal
        x = zeros(n,1);
        T = randi([1,k]);
        q = randperm(T);
        x(q(1:T)) = sign(randn(T,1)); 

        % Calculando a variaância do ruído de acordo com o SNR
        sigma2 = norm(x)*10^(-SNR/10);
        
        % Vetor ruído
        noise = sqrt(sigma2)*randn(m, 1);
        
        % observations
        y = A*x + noise;
        
        % initial guess = min energy
        x0 = A'*y;
        
        % solve the LP
        tic
        xp = l1eq_pd(x0, A, [], y, 1e-4);
        toc
        
        %Medida de sucesso
        erro = norm(x-xp);
        if(erro <= Bound_S2)
            acerto = acerto + 1;
        end
    end
    
    if (acerto/Iterations > Bound_S1)
        k = k + 1;
        acerto_anterior = acerto/Iterations;
    else
        SucessMap(m, k) = acerto_anterior;
        %keyboard()
        k = 1;
        m = m + HorizontalDistance;   
    end
    %keyboard()
end

eixo_x = 1/n:1/n:1;
eixo_y = 1:-1/n:1/n;

%save('SucessNoNoise.mat','SucessMap');
save('SucessMapSNR75.mat','SucessMap');

heatmap(eixo_x, eixo_y, rot90(SucessMap))