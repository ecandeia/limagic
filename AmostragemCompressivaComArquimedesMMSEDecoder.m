%%
% Evaluate the mean norm evoluation.
% Gaussin Source N(0,1)
% AWGN Channel according CSNR (dB)
% Receiver MMSE, quantize the output, table search to determine input
% that produces nearest output to wuantized value
% Input:
%   Delta: radial distance betwwen to neighboor spiral arms (it needs
%          to be optimized.
% Output:
%   Vector with the evoluation of the norm.
%   Mean norm for a given CSNR
%
% Edmar C. Gurjao -- ecandeia@dee.ufcg.edu.br 
% 03/22/2012
%%

clear all;

% put key subdirectories in path if not already there
path(path, './Optimization');
path(path, './Data');

% Parâmetros:

    %Número de pontos na curva
    pts_curva = 45;
    
    %Contador
    contador=75;

    %Delta = [2.15 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05];
    Delta = ones(1,pts_curva+1)/20;
    % Table of quantized values to be used in the reception
    x_temp = -5:0.0001:5;

    % Channel signal to Noise Ratio (dB)
    pos = 1;
    count=0;
    sigma_s1 = 0;

for CSNR = 0:1:pts_curva 
    
    %Jow(CSNR+1)= Delta (pos);
    %Jon(pos) = Delta(pos);
    
    %Curva de arquimedes table(1,:),table(2,:)
    table = (1/pi).*[sign(x_temp).*Delta(pos).*sqrt(6.25*abs(x_temp)/Delta(pos)).*cos(sqrt(6.25*abs(x_temp)/Delta(pos)));
    sign(x_temp).*Delta(pos).*sqrt(6.25*abs(x_temp)/Delta(pos)).*sin(sqrt(6.25*abs(x_temp)/Delta(pos)))];
    
    %table = (1/pi).*[sign(x_temp).*Delta.*sqrt(6.25*abs(x_temp)/Delta).*cos(sqrt(6.25*abs(x_temp)/Delta));
    %sign(x_temp).*Delta.*sqrt(6.25*abs(x_temp)/Delta).*sin(sqrt(6.25*abs(x_temp)/Delta))];

    M = 2;
    N = 20;
    S = 2;
    
    % measurement matrix
    Phi = randn(M,N);
    Phi = orth(Phi')';
    
    
    for repeat = 1:contador
        % Gaussian source N(0,1) -> Fonte Gaussiana com média 0 e variância 1
        x = zeros(N,1);
        q = randperm(N);
        x(q(1:S)) = 0.25*randn(S,1);
     
        count=count+1;
        
        %ponto transmitido
 %       x_fonte(count) = x;
                
        y = Phi*x;

        % Ponto a ser transmitido pela fonte
        z(1,:) = (1/pi).*[sign(y(1,1)).*Delta(pos).*sqrt(6.25*abs(y(1,1))/Delta(pos)).*cos(sqrt(6.25*abs(y(1,1))/Delta(pos)));
    sign(y(1,1)).*Delta(pos).*sqrt(6.25*abs(y(1,1))/Delta(pos)).*sin(sqrt(6.25*abs(y(1,1))/Delta(pos)))];

        z(2,:) = (1/pi).*[sign(y(2,1)).*Delta(pos).*sqrt(6.25*abs(y(2,1))/Delta(pos)).*cos(sqrt(6.25*abs(y(2,1))/Delta(pos)));
    sign(y(2,1)).*Delta(pos).*sqrt(6.25*abs(y(2,1))/Delta(pos)).*sin(sqrt(6.25*abs(y(2,1))/Delta(pos)))];
        
        %Valor enviado pela fonte, tá na curva de arquimedes
%         z1c (count) = z(1);
%         z2c (count) = z(2);
        
        
        % Var é  a variancia
        sigma_z = var(z)';
        %sigma_z1 (count) = var(z);
        
        % Acrescenta o Ruido (n)
%         ruido = randn(2,1);
%         ruido1c(count) = ruido(1);
%         ruido2c(count) = ruido(2);
%         n = sqrt(sigma_z*10^(-CSNR/10)/2)*ruido;
%         n1c (count) = n(1);
%         n2c (count) = n(2);
        n(1,:) = sqrt(sigma_z(1,1)*10^(-CSNR/10)/2)*randn(2,1);
        n(2,:) = sqrt(sigma_z(2,1)*10^(-CSNR/10)/2)*randn(2,1);
        
        %R é a soma sinal + ruido, Sinal enviado + ruido 
        zhat = z;% + n; %sinal sem quantização
%         r1c(count) = zhat(1) + n(1); 
%         r2c(count) = zhat(2) + n(2);
        
        %Faz quantização do sinal enviado+ruido
        [index_1,zhat_quant(1,:)] = quantiz(zhat(1,:),-0.900:0.001:0.200,-0.900:0.001:0.201); %sinal quantizado
        [index_2,zhat_quant(2,:)] = quantiz(zhat(2,:),-0.900:0.001:0.200,-0.900:0.001:0.201); %sinal quantizado
        
%         zhat_quant1c (count)= zhat_quant(1);
%         zhat_quant2c (count) = zhat_quant(2);
%         index1c (count) = index(1);
%         index2c (count) = index(2);
        
        indice_1 = 1;
        colnorm_1 = sqrt(sum((table - repmat(zhat_quant(1,:)',1,size(table,2))).^2,1)); % Distâcia entre os pontos recebidos

        indice_2 = 1;
        colnorm_2 = sqrt(sum((table - repmat(zhat_quant(2,:)',1,size(table,2))).^2,1)); % Distâcia entre os pontos recebidos
        
        [r_tmp(1,1) indice_1] = min(colnorm_1); %Menor valor da distância
        [r_tmp(2,1) indice_2] = min(colnorm_2); %Menor valor da distância
        %Sinal Estimado, pela menor distância (MMSE método)
        
        yhat(1,1) = x_temp(indice_1);
        yhat(2,1) = x_temp(indice_2);
        
%         sinal_estimado1c(count) = table(1,indice);
%         sinal_estimado2c(count) = table(2,indice);

        xhat0 = Phi'*yhat;
        
        % solve the LP
        xhat = l1eq_pd(xhat0, Phi, [], yhat, 1e-3);
        
        %O que faz daqui pra baixo?
        norm_evol_MMSE = norm(x-xhat);
        evol_MMSE_1c(count) = norm_evol_MMSE;
    end
    var_SDR(pos) = var(norm_evol_MMSE) + mean(norm_evol_MMSE)^2;
    %até aqui, não entendi....
    
    pos = pos+1;
end
% s_novo = [s1c;s2c]; %sinal transmitido
% n_novo = [n1c;n2c]; %ruido adquirido
% r_novo = [r1c;r2c]; %sinal (ruido + sinal)
% r_quantizado = [r_quant1c;r_quant2c]; %sinal quantizado

% stem(x);
% hold on 
% stem(xhat,'r')
%%

CSNR = 0:1:pts_curva;
CSNR2 = 10.^(0.1*CSNR); %Para calcular o OPTA 1:2
SNR_theo = (1+CSNR2).^2; %Para calcular o OPTA 1:2

SDR_MMSE = 10*log10(0.25^2./var_SDR); % O que é essa função?

%começo da amostragem compressiva
% A = randn(pts_curva+1,1);
% y = A*SDR_MMSE;
% x0 = A.'*y;
%xp = l1eq_pd(x0, A, [], y);
%fim da amostragem (falta o calculo de l1eq),
figure (1);
grid on;
hold on;
graph1 = plot(CSNR,SDR_MMSE, 'r');
graph2= plot(CSNR,10*log10(SNR_theo),'k'); %Gráfico OPTA 1:2;
%graph3 = plot (xp,'b'); %Grafico da Amostragem compressiva de SDR_MMSE
legend([graph1 graph2],'SDR_MMSE - Curva do MMSE', 'Gráfico OPTA 1:2');
%legend([graph1 graph2 graph3],'SDR_MMSE - Curva do MMSE', 'Gráfico OPTA 1:2','Xp - Amostragem Compressiva de SDR_MMSE');
xlabel('CSNR (dB)');
ylabel('SDR_MMSE');

% figure (2);
% grid on;
% hold on;
% grafico1= plot (table(1,:),table(2,:),'g'); %Curva de Arquimedes (Y)
% grafico2= plot (s1c,s2c,'x'); %Pontos que a fonte envia (Z)
% grafico3 =plot (r1c,r2c,'*c'); %Sinal recebido sem quantização 
% grafico4= plot (r_quant1c,r_quant2c,'o'); %Sinal Recebido Quantizado 
% grafico5 =plot (sinal_estimado1c,sinal_estimado2c,'d'); %Sinal Estimado (Y chapeu)
% legend([grafico1 grafico2 grafico3, grafico4 grafico5],'Curva de Arquimedes', 'Amostras da fonte', 'Sinal Recebido sem quantização', 'Sinal Recebido Quantizado','Sinal Estimado')
% axis([-0.3 0.3 -0.3 0.3])
%plot (s1c,s2c,'c*');

%save spiral_MMSE SDR_MMSE;