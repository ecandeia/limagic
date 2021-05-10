% l1dantzig_example.m
%
% Test out l1dantzig code (l1 minimization with bounded residual correlation).
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

% Faz a busca pelo valor de epsilon (restrição da otimização que fornece o
% menor erro. No momento faz a busca exaustiva, partindo de um ponto
% aleatoriamente escolhido, até atingir o menor valor de erro. Nesse ponto
% para e guarda os valores.
%
% Uma possibilidade é usar a técnica descrita em: PRIMAL DUAL PURSUIT A 
% HOMOTOPY BASED ALGORITHM FOR THE DANTZIG SELECTOR, Muhammad Salman Asif.
% Georgia Institue of Technology, 2008.

clear all

% put optimization code in path if not already there
path(path, './Optimization');

% signal length
N = 80;

% number of spikes to put down
T = 5;

% number of observations to make
K = 30;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1));

% measurement matrix = random projection
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
disp('Done.');
	
% noisy observations
sigma = 0.5;
e = sigma*randn(K,1);
y = A*x + e;

% initial guess = min energy
x0 = A'*y;

% Dantzig selection

lambda = 0.01:0.01:2;
erro_ant = Inf;
for index = 1:length(lambda)
    xp = l1dantzig_pd(x0, A, [], y, lambda(index), 5e-2);
    erro(index) = norm(x-xp);
    if erro(index) < erro_ant
        eps = lambda(index);
        erro_ant = erro(index);
    else
        break;
   end
end



% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1dantzig_pd(x0, Afun, Atfun, y, epsilon, 5e-2, 50, 1e-8, 500);
% toc

plot(lambda(1:index), erro(1:index))



