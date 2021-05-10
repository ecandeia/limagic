%!= Using Compressed Sensing to transmit in a AWGN Channel with spike noise.
%   Recovey by norm minization via CVX
%
% Last verstion 09/07/2012
%
% Edmar Candeia Gurjão
% Altered by Edmar at 14:00h (USA time)
%
%=!

clear all
tic

TSNR_MAX  = 30;
TSNR_MIN  = 0;
TSNR_INCR = 5;

TSNR_dB = TSNR_MIN:TSNR_INCR:TSNR_MAX;

D_est = zeros(length(TSNR_dB),1);
    
% Source variance
sigma_x2 = 1;


% Parameters for Compressive Sensing

% signal length
N = 100;

% number of observations to make
M = 80;

% Spartsity of the source vector
T = 6;

% Sparsity of the spike noise
S = 0;

% Number os samples vectors
Nrepetitions = 1e2;


% Total Power
Pt = 1;

for k = 1:length(TSNR_dB)
    
    str = sprintf('CSNR %d (dB)', TSNR_dB(k));
    disp(str);
    
    TSNR = 10^(0.1*TSNR_dB(k));  

    sigma_n2 = Pt/(TSNR);    
    
    %A=sqrt(1/M)*randn(M,N);
    %W=diag(1./sqrt(diag(A'*A)));
    %A=A*W;
    A = randn(M,N);
    A = orth(A')';
      
    i4 = 1;
    i5 = 1;
    for nr = 1:Nrepetitions
        
        % Sparse Source Symbols
        x = zeros(N,1);
        q = randperm(N);
        x(q(1:T)) = sqrt(sigma_x2)*randn(T,1);

        nx(i5) = norm(x)^2;
        i5 = i5 +1;
        
        % Compress
        xc = A*x;
           
        gamma = sqrt(Pt/(var(xc) + mean(xc)^2));
        
        xc = gamma*xc;

        %% Noise
        % Gaussian noise
        noise = sqrt(sigma_n2)*(randn(size(xc,1),1));     
        
        beta = gamma*(var(xc) + mean(xc)^2)/(gamma^2*(var(xc) + mean(xc)^2) + sigma_n2);
        
        % Spike noise
        s = zeros(M,1);
        qs = randperm(S);
        s(qs(1:S)) = sign(randn(S,1));
        
        
        zhat = xc + noise + s; 
                
        zhat = beta*zhat;
        
        %% Recovery from compressed measurements
%         cvx_begin
%            variable xhat(N);
%            minimize(norm(A*xhat-zhat,1));
%         cvx_end;
% Dantzig selection
epsilon = 3e-3;
x0 = A'*zhat;
xhat = l1dantzig_pd(x0, A, [], zhat, epsilon, 5e-2);

%         support = find(abs(xhat)>0.2);
%         
%         AS=A(:,support);
%         Qs=AS'*AS/sigma_n2^2;% + eye(T)/sigma_x2;
%         IQs=inv(Qs);
%         x_oracle=zeros(N,1);
%         
%         x_oracle(support)=(Qs\AS')*zhat/sigma_n2^2;%IQs*AS'*zhat/sigma_n2^2;
%         
        xhat_final = xhat;%x_oracle;

        %% Estimated Distortion  
  
        for i3 = 1:N
         pontual_distortion(i4,1) = (norm(x(i3,:) - xhat_final(i3,:))^2);
         i4 = i4 +1;
        end
        
    end
    
   D_est(k) = mean(pontual_distortion);
    
 
end



%Plot %CSNR x Distortion 

%CSNR = 10.^(0.1*CSNR_dB);
%TSNR = CSNR*M;
%TSNR_dB = 10*log10(TSNR);

if (M == 20)
    plot(TSNR_dB,10*log10(mean(nx)./D_est),'k-*');
else
    plot(TSNR_dB,10*log10(mean(nx)./D_est),'r-*');
end
    
xlabel('TSNR');
ylabel('SDR');
grid on
toc

%legend('OPTA','Calculated','Simulation');