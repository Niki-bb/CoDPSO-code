function [swarm]=inicializaSwarm(N, N_PAR, V_MAX, V_MIN)


% warning off last      %impede a escrita do ultimo warning dado (variaveis globais)
% global N;
% global N_PAR;
% global V_MIN;
% global V_MAX;

%% ÐÞ¸ÄÇ°
  swarm = zeros(N,N_PAR);
  for i = 1: N             %varia de 1 at? ao n? de testes pretendido(N) 
      for j = 1: N_PAR     %dentro de cada teste vai variar de 1 at? 4(n? de parametros k)
          swarm(i,j) = rand(1,1) * ( V_MAX-V_MIN ) + V_MIN;   %cria uma matriz ixj com os valores de 
      end                                                        %todos os k
  end