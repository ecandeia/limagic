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
clear all

% signal length
N = 1000;

% number of observations to make
M = 400;

% Encoding size
K = 600;

% number of spikes in the signal
T = 20;

% Parity size
km = 200;

% Errors
TE = 10;

% random +/- 1 signal
x = zeros(N,1);
q = randperm(N);
x(q(1:T)) = sign(randn(T,1));

% measurement matrix
disp('Creating measurment matrix...');
D = randn(M,N);
D = orth(D')';
disp('Done.');

% observations
y = D*x;

% encoding matrix
disp('Creating encoding matrix...');
H = randn(km,K)/sqrt(M);
disp('Done.');

tmp = null(H);

G = tmp(:,1:M);

% encoding
c = G*y;

% Errors
e = zeros(K,1);
q = randperm(M);
e(q(1:TE)) = randn(TE,1);


r = c + e;

s = H*r;

% recover
r0 = H'*s;%inv(H'*H)*H'*s;
ehat = l1decode_pd(r0, H, [], s, 1e-4, 30);


c_hat = c - ehat;

% initial guess = min energy
c_hat_0 = D'*c_hat;

% solve the LP
tic
xp = l1eq_pd(c_hat_0, D, [], c_hat, 1e-3);
toc

% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);
% toc




