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

% To reproduce the example in the documentation, uncomment the 
% two lines below
%load RandomStates
%rand('state', rand_state);
%randn('state', randn_state);

% signal length
N = 500;
% number of spikes in the signal
T = randi([1,5]);
% number of observations to make
K = 20;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(T);
x(q(1:T)) = sign(randn(T,1));


% measurement matrix
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
disp('Done.');

% observations
y = A*x;


%stem(y)

% initial guess = min energy
x0 = A'*y;

% solve the LP
tic
xp = l1eq_pd(x0, A, [], y, 1e-3);
toc

%figure(2)
stem(x)
hold on

%figure(2)
stem(xp,'r')

norm(x-xp)

% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);
% toc




