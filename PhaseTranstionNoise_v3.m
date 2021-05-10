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

ProblemSets = 2;
HorizontalDistance = 2;
Bound_S1 = 0.9;%1e-3;
Bound_S2 = 1e-3;

k = 1;

n = 100;

m = HorizontalDistance;

CountMap = zeros(n, n);

SucessMap = zeros(n,n);

Iterations = 100;
acerto_anterior = 0;
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
        
        % observations
        y = A*x;
        
        % initial guess = min energy
        x0 = A'*y;
        
        % solve the LP
        tic
        xp = l1eq_pd(x0, A, [], y, 1e-4);
        toc
        
        %Medida de sucesso
        erro = norm(x-xp);
        if(erro <= 1e-2)
            acerto = acerto + 1;
        end
    end
    
    if (acerto/Iterations > 0.9)
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

heatmap(eixo_x, eixo_y, rot90(SucessMap))