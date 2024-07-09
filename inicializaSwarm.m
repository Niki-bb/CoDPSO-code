function [swarm]=inicializaSwarm(N, N_PAR, V_MAX, V_MIN)


% warning off last      %impede a escrita do ultimo warning dado (variaveis globais)
% global N;
% global N_PAR;
% global V_MIN;
% global V_MAX;

%% –ﬁ∏ƒ«∞
  swarm = zeros(N,N_PAR);
  for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
      for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
          swarm(i,j) = rand(1,1) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
      end                                                        %todos os k
  end
end
%   
% 
%  %% 1.1 ChebyshevªÏ„Á
%  Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=cos(count*acos(Z(count)));
%           end
%       end                                                        %todos os k
%   end
  
%  %% 1.2 CircleªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=mod(Z(count)+0.2-(0.25/pi)*sin(2*pi*Z(count)),1);
%           end
%       end                                                        %todos os k
%   end  
  
% %    %% 1.3 Gauss/mouseªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               if Z(count)==0
%                   Z(count+1)=0;
%               else
%                  Z(count+1)=1/mod(Z(count),1);
%               end
%           end
%       end                                                        %todos os k
%   end  
  
%     %% 1.4  IterativeªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=sin(1*pi/Z(count));
%           end
%       end                                                        %todos os k
%   end  
 
  
% %% 1.5 ¬ﬂº≠ªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=0.7;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=2*Z(count)*(1-Z(count));
%           end
%       end                                                        %todos os k
%   end

%     %% 1.6  Piecewise linearªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% P=0.4;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               if 0 <= Z(count) < P
%                     Z(count+1)=Z(count)/P;
%               elseif P <= Z(count) < 0.5
%                     Z(count+1)=(Z(count)-P)/(0.5-P);
%               elseif 0.5 <= Z(count) < (1-P)
%                      Z(count+1)=(1-Z(count)-P)/(0.5-P);
%               elseif (1-P) <= Z(count) < 1
%                     Z(count+1)=(1-Z(count))/P;
%               end
%           end
%       end                                                        %todos os k
%   end  


%     %% 1.7  SineªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% A8=4;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=(A8/4)*sin(pi*Z(count));
%           end
%       end                                                        %todos os k
%   end  
%   
 
%     %% 1.8  SingerªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% miu=1.07;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=miu*(7.86*Z(count)-23.31*Z(count)^2+28.75*Z(count)^3-13.3*Z(count)^4);
%           end
%       end                                                        %todos os k
%   end  
  
%      %% 1.9  SinusoidalªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=2.3*(Z(count)^2)*sin(pi*Z(count));
%           end
%       end                                                        %todos os k
%   end  

%    %% 2.0 TentªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% b6=0.5;
% b8=2;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               if 0 <= Z(count) < 0.5
%                   Z(count+1)=Z(count)/b6;
%               elseif 0.5 <= Z(count) <= 1
%                  Z(count+1)=b8*(1-Z(count));
%               end
%           end
%       end                                                        %todos os k
%   end  


%        %% 2.1-1  LCSªÏ∫œªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% r=0.7;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               Z(count+1)=cos(pi*(4*r*Z(count)*(1-Z(count))+(1-r)*sin(pi*Z(count))-0.5));
%           end
%       end                                                        %todos os k
%   end  

%      %% 2.1-2  LCSªÏ∫œªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% r=0.7;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               if Z(count) < 0.5
%                 Z(count+1)=cos(pi*(r*sin(pi*Z(count))+2*(1-r)*Z(count)-0.5));
%               else
%                 Z(count+1)=cos(pi*(r*sin(pi*Z(count))+2*(1-r)*(1-Z(count))-0.5)); 
%               end
%           end
%       end                                                        %todos os k
%   end 

%      %% 2.1-3  LCSªÏ∫œªÏ„Á
% Z=zeros(1,N*N_PAR);
% Z(1)=rand;
% r=0.7;
% swarm = zeros(N,N_PAR);
%   for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
%       for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
%           for count=1:N*N_PAR
%               swarm(i,j) = Z(count) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
%               if Z(count) < 0.5
%                 Z(count+1)=cos(pi*(2*r*Z(count)+4*(1-r)*Z(count)*(1-Z(count))-0.5));
%               else
%                 Z(count+1)=cos(pi*(2*r*(1-Z(count))+4*(1-r)*Z(count)*(1-Z(count))-0.5)); 
%               end
%           end
%       end                                                        %todos os k
%   end 
% 

%% –ﬁ∏ƒ«∞