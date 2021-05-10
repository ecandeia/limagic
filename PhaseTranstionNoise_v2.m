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
HorizontalDistance = 10;
Bound_S1 = 0.9;%1e-3;
Bound_S2 = 1e-6;

k = 1;

n = 100;

m = HorizontalDistance;

CountMap = zeros(n, n);

SucessMap = zeros(n,n);

%TL = zeros(n);

while (m < n)

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
    xp = l1eq_pd(x0, A, [], y, 1e-3);
    toc
                
    CountMap(k, m) = CountMap(k, m) + 1;
    
    %Medida de sucesso
    CheckOnSucess = norm(x-xp)/norm(x);
    
    if (CheckOnSucess < Bound_S2)
        SucessMap(k, m) = SucessMap(k, m) + 1;
        disp("Sucesso")
        keyboard()
    end
    %keyboard()
    if ((CountMap(k, m) >= ProblemSets)&&((SucessMap(k, m)/CountMap(k, m) >= Bound_S1)||k==1))
        if ((k == 1)&&(SucessMap(k, m)/CountMap(k, m) >= Bound_S1))
            TL(m) = 0;
        else
            TL(m) = k/m;
        end
        %keyboard()
        disp("If 1")
        m = m + HorizontalDistance;
    elseif (SucessMap(k, m)/CountMap(k, m) >= Bound_S1)
        if (k < m) 
            k = k + 1;
        elseif ((SucessMap(k, m)/CountMap(k, m) < Bound_S1)&&( k > 1))
            k = k - 1;
        end
    end
    %keyboard()
end

indice = 1;
for i = HorizontalDistance:HorizontalDistance:n - HorizontalDistance
    k = round(i*TL(i))
    grafico(k, i) = TL(i);
end
   
stem(grafico)