%%%%%%%%%%%%%%%%%2024.2.12 wxb修改注释 邮箱：649713453@qq.com%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;              %清除所有变量
% tic                     %记录程序运行时间  计时开始
% clc;                    %清屏
function [Convergence_curve,data]=CoDPSO(fhd,pop_size1,dim,Max_iteration,lb,ub,aer,bei,varargin)
fbias=[300, 400, 600, 800, 900, 1800,...
       2000, 2200, 2300, 2400, 2600, 2700];
population=pop_size1;         %群体粒子个数
N_PAR=dim;               %粒子维数
N_GER=Max_iteration;               %最大迭代次数
nswarms=5;              %群体数
W=1;
PHI1=1.5;               %学习因子1
PHI2=1.5;               %学习因子2
% [Us]=textread('25DST电压.txt','%.2f');    
% Ut=(Us(2430:12230))';
aerfa=aer;
beita=bei;
Convergence_curve=zeros(1,N_GER); 
data = [];

N_MIN=10;            % 10 本来是round(population/2)
N_Init=population;                    % 30
N_MAX=50;                             % 本来是population*2
MIN_SWARMS=round(nswarms/2);          % round函数表示四舍五入取整 3
N_SWARMS=nswarms;                     % 5
MAX_SWARMS=nswarms*2;                 % 10

SWARMS=zeros(MAX_SWARMS,1);           % (10,1)
for i=1:MAX_SWARMS
    if (i<=N_SWARMS)
        SWARMS(i)=1;                  % SWARMS=[1;1;1;1;1;0;0;0;0;0]
    end
end

N=N_Init*ones(MAX_SWARMS,1);          % 维度为[10,1],其中元素值均为N_Init-20
% SCmax = round(N_GER/10);               % 关于SCmax的取值多参考一些文献，自己调试一下 此处SCmax为1000
SCmax = 30; 
SC=zeros(MAX_SWARMS,1);               % (10,1)

X_MAX=ub;
X_MIN=lb;
vmax = (X_MAX-X_MIN)/4;
vmin = -vmax;

gbestvalue = 1000000*ones(MAX_SWARMS,1); 
gbestvalueaux=gbestvalue;

fit = cell(MAX_SWARMS,1);       % 均为(10,1) cell数组一般被叫做元胞数组，它的每个单元可以储存不同的数据类型，可以是数值，字符或矩阵或元胞数组等，类似于学过的c语言里的结构体。本程序中元胞数组的每个单元储存的是矩阵
x=cell(MAX_SWARMS,1);
v=cell(MAX_SWARMS,1);
vbef1=cell(MAX_SWARMS,1);
vbef2=cell(MAX_SWARMS,1);
vbef3=cell(MAX_SWARMS,1);
for i=1:N_SWARMS                % i=1:5
    x{i}=zeros(N(i),N_PAR);     % fitness do vector best (20,N_PAR)全零阵   
    fit{i}(1:N(i),1)=0;         % fitness do vector best 
    v{i}=zeros(N(i),N_PAR);     % inicializa v a NxN_PAR zeros
    vbef1{i}=zeros(N(i),N_PAR);  % 具体问题应该初始化为不同的范围
    vbef2{i}=zeros(N(i),N_PAR);
    vbef3{i}=zeros(N(i),N_PAR);
end
xaux=x;
fitBest=fit;
vaux=v;
vbef1aux=vbef1;       % 储存历代速度值已进行分数阶速度更新
vbef2aux=vbef2;
vbef3aux=vbef3;

xBest=cell(MAX_SWARMS,1);
gBest=cell(MAX_SWARMS,1);
Nkill=zeros(MAX_SWARMS,1);   %(10,1)全零阵
nger=1;
%%%%%%%%%%%%%%%%%%初始化个体最优位置和最优值%%%%%%%%%%%%%%%%%%%
for i=1:N_SWARMS        % i=1:5
    xaux{i}(1:N(i),1:N_PAR)=inicializaSwarm(N(i), N_PAR, X_MAX, X_MIN); %x ser? uma matriz com todos os valores de k  / inicializaSwarm译为初始化粒子群，不太懂这一步含义是啥？
    x=xaux;
    for j=1:N(i)
        fit{i}(j)=feval(fhd,(x{i}(j,:))',varargin{:})-fbias(varargin{:});
        fitBest{i}(j)=fit{i}(j);      %INICIALIZAR FITBEST  fitBest储存每个群体中每个粒子的适应度值
    end
    
    [a,b]=min(fit{i});
    gBest{i}=x{i}(b,:);           % gBest储存每个群体中的最优粒子值
    gbestvalue(i) = fit{i}(b);    % gbestvalue储存每个群体中最优粒子的适应度值
    gbestvalueaux(i)=gbestvalue(i);  %
    % Aqui inicializa-se o xBest (Matriz do melhor local para cada part韈ula)
    xBest{i}(1:N(i),1:N_PAR) = inicializaSwarm(N(i), N_PAR, X_MAX, X_MIN); %melhor x
end

fitaux=fit;
xBestaux=xBest;
gBestaux=gBest;
fitBestaux=fitBest;  
gb=ones(N_GER,1);                % 记录迭代过程中每一代的最优适应度值
while (nger<=N_GER)
    for i=1:MAX_SWARMS      % i=1:10
        if (SWARMS(i)==1)   %swarm is alive
            
            xBestaux{i}=xBest{i};
            vbef1aux{i}=vbef1{i};
            vbef2aux{i}=vbef2{i};
            vbef3aux{i}=vbef3{i};
            randnum1 = rand ([N(i), N_PAR]);
            randnum2 = rand ([N(i), N_PAR]);
            
            % 速度更新
            gaux=ones(N(i),1);
            vaux{i} = W*(aerfa*v{i}(1:N(i),:) - (1/2)*(aerfa^2-aerfa-beita^2)*vbef1{i}(1:N(i),:)...
                      + (1/6)*(aerfa^3-3*aerfa^2+aerfa*(2-3*beita^2)+3*beita^2)*vbef2{i}(1:N(i),:)...
                      + (1/24)*(aerfa^4+4*beita^4-4*aerfa^3-6*aerfa^2*beita^2+11*(aerfa^2-beita^2)+18*aerfa*beita^2-6*aerfa)*vbef3{i}(1:N(i),:))...
                      + randnum1.*(PHI1.*(xBest{i}(1:N(i),:)-x{i}(1:N(i),:)))...
                      + randnum2.*(PHI2.*(gaux*gBest{i}-x{i}(1:N(i),:)));
            % 速度限制
            for j=1:N_PAR
                for q=1:N(i)
                    vaux{i}(q,j) = ( (vaux{i}(q,j) <= vmin).*vmin ) + ( (vaux{i}(q,j) > vmin).*vaux{i}(q,j) );
                    vaux{i}(q,j) = ( (vaux{i}(q,j) >= vmax).*vmax ) + ( (vaux{i}(q,j) < vmax).*vaux{i}(q,j) );
                end
            end
            vbef3aux{i}=vbef2{i};
            vbef2aux{i}=vbef1{i};
            vbef1aux{i}=v{i};
            % 位置更新
            xaux{i} = x{i}(1:N(i),:)+vaux{i}(1:N(i),:);

            % limit the position between the mimimum and maximum thresholds
            % 位置限制
            for j = 1:N(i)
                for k = 1:N_PAR
                    if xaux{i}(j,k) < X_MIN
                        xaux{i}(j,k) = X_MIN;
                    elseif xaux{i}(j,k) > X_MAX
                        xaux{i}(j,k) = X_MAX;
                    end
                end
            end
            
            for j=1:N(i)
                fitaux{i}(j)=feval(fhd,(xaux{i}(j,:))',varargin{:})-fbias(varargin{:});
                if fitaux{i}(j) < fitBestaux{i}(j)     % fitBestaux是存储的每个群体中每个粒子的适应度值，相当于初始化fitaux的值，fitaux是当前的。
                   fitBestaux{i}(j) = fitaux{i}(j);
                   xBestaux{i}(j,:) = xaux{i}(j,:);
                end
            end
            
            [a,b]=min(fitaux{i});                 % 或许这里应该是min(fitBestaux)?
            if (a < gbestvalueaux(i))
                gBestaux{i}=xaux{i}(b,:);         % gBestaux储存每个群体中的最优粒子
                gbestvalueaux(i) = fitaux{i}(b);  % gbestvalueaux储存每个群体中的最优粒子的适应度值
            end
            [fitting(nger),indexM]=min(gbestvalueaux);  % fitting储存所有群体中的最优适应度值
%             gBestT(nger,1:N_PAR) = gBestaux{indexM};    % gBestT储存所有群体中具有最优适应度值的粒子
            
            if (nger>1)
                if gbestvalueaux(i)>=gbestvalue(i)  %check performance of the swarm  这里应该是考虑最小值问题时，最大值问题这里应该是小于等于，后续程序进行自适应调整
%                     if gbestvalueaux(i)>2870      %
%                     DE策略，2870为退化点，根据算法表现判断是否需要添加此策略
                        SC(i)=SC(i)+1;
                        if (SC(i)==SCmax)        %reached search limit  SCmax
                            if N(i)>N_MIN        %delete particle  N_MIN
                               [a,b]=max(fitaux{i});    % 因为每次惩罚值删除一个粒子，所以这里只考虑一个b（因为可能存在多个粒子都具有最坏的适应度值，只是目标解的维度较大时这个可能比较小）

                                if b==1      % 当这个最大值在第一行时（假设是最小值问题）
                                    xaux{i}=xaux{i}(2:N(i),:);
                                    xBestaux{i}=xBestaux{i}(2:N(i),:);
                                    vaux{i}=vaux{i}(2:N(i),:);
                                    vbef1aux{i}=vbef1aux{i}(2:N(i),:);
                                    vbef2aux{i}=vbef2aux{i}(2:N(i),:);
                                    vbef3aux{i}=vbef3aux{i}(2:N(i),:);
                                end
                                if b==N(i)   % 当这个最大值在最后一行时
                                    xaux{i}=xaux{i}(1:b-1,:);
                                    xBestaux{i}=xBestaux{i}(1:b-1,:);
                                    vaux{i}=vaux{i}(1:b-1,:);
                                    vbef1aux{i}=vbef1aux{i}(1:b-1,:);
                                    vbef2aux{i}=vbef2aux{i}(1:b-1,:);
                                    vbef3aux{i}=vbef3aux{i}(1:b-1,:);
                                end
                                if (b>1 && b<N(i))
                                    xaux{i}=vertcat(xaux{i}(1:(b-1),:),xaux{i}((b+1):N(i),:));    % vercat表示垂直连接数组
                                    xBestaux{i}=vertcat(xBestaux{i}(1:(b-1),:),xBestaux{i}((b+1):N(i),:));
                                    vaux{i}=vertcat(vaux{i}(1:(b-1),:),vaux{i}((b+1):N(i),:));
                                    vbef1aux{i}=vertcat(vbef1aux{i}(1:(b-1),:),vbef1aux{i}((b+1):N(i),:));
                                    vbef2aux{i}=vertcat(vbef2aux{i}(1:(b-1),:),vbef2aux{i}((b+1):N(i),:));
                                    vbef3aux{i}=vertcat(vbef3aux{i}(1:(b-1),:),vbef3aux{i}((b+1):N(i),:));
                                end
                                N(i)=N(i)-1;
                                Nkill(i)=Nkill(i)+1;
                                SC(i)=fix(SCmax*(1-1/(Nkill(i)+1)));     % fix函数表示向零方向取整
                            else                    %delete swarm
                                if (N_SWARMS>MIN_SWARMS)
                                    SWARMS(i)=0;
                                    N_SWARMS=N_SWARMS-1;
                                    SC(i)=0;
                                    N(i)=0;
                                    %                                         fprintf('\ndelete swarm %d',i);
                                end
                            end
                        end
%                     end
                else
                    if (Nkill(i)>0)           % Nkill初始值为0
                        Nkill(i)=Nkill(i)-1;
                    end
                    if (N(i)<N_MAX)     %create particle  N_MAX=40
                        
                        N(i)=N(i)+1;
                        %N_ALIVE=size(xaux{i},1); 
                        xaux{i}(N(i),1:N_PAR)=inicializaSwarm(1, N_PAR, X_MAX, X_MIN);     % 为新加粒子进行赋值
                        xBestaux{i}(N(i),1:N_PAR)=inicializaSwarm(1, N_PAR, X_MAX, X_MIN); % 初始化新加粒子的最优位置
                        vaux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);     % 初始化新加粒子的速度值
                        vbef1aux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);
                        vbef2aux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);
                        vbef3aux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);
                        fitaux{i}(N(i))=feval(fhd,(xaux{i}(N(i),:))',varargin{:})-fbias(varargin{:});
                        fitBestaux{i}(N(i))=fitaux{i}(N(i)); % 储存新加粒子的适应度值
                        [a,b]=min(fitaux{i});     % 找出新加粒子后的种群中的最小适应度值，此时新加粒子的适应度值已经包含在fitaux中
                        gBestaux{i}=xaux{i}(b,:); % 储存新加粒子后的种群的最优位置
                        gbestvalueaux(i) = fitaux{i}(b); % 储存新加粒子后的种群的最优适应度值
                        %                                 fprintf('\ncreate particle in swarm %d',i);
                    end
                end
                if (Nkill(i)==0)&&(N_SWARMS<MAX_SWARMS)     %add new swarm N_SWARMS初始值为5，MAX初始值为10
                    
                    prob=rand/N_SWARMS;
                    create_swarm=rand;
                    if (create_swarm<prob)
                        SC(i)=0;                        % SC初始值为0
                        SWARM_ALIVE=find(SWARMS==0);    % search the swarms that are dead 此行表示返回SWARMS中等于0的元素的位置
                        SWARMS(SWARM_ALIVE(1))=1;       % 将第一个等于0的SWARMS赋值为1（原来为0）
                        N_SWARMS=N_SWARMS+1;
                        N(SWARM_ALIVE(1))=N_Init;       % N_Init为种群大小，等于20
                        
                        %%% 修改补充
                        [~,b]=sort(fitaux{i});
                        xaux{SWARM_ALIVE(1)}(1:N_Init/2,:)=xaux{i}(b(1:N_Init/2),:);  % 初始化新加种群的位置
                        vaux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vaux{i}(b(1:N_Init/2),:);       % 初始化新加种群的速度值
                        vbef1aux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vbef1aux{i}(b(1:N_Init/2),:);
                        vbef2aux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vbef2aux{i}(b(1:N_Init/2),:);
                        vbef3aux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vbef3aux{i}(b(1:N_Init/2),:);
                        xBestaux{SWARM_ALIVE(1)}(1:N_Init/2,:) = xBestaux{i}(b(1:N_Init/2),:);  % 初始化新加种群的历史最优位置
                        %%% 修改补充
                        
                        xaux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=inicializaSwarm(N_Init/2, N_PAR, X_MAX, X_MIN);  % 初始化新加种群的位置
                        vaux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);       % 初始化新加种群的速度值
                        vbef1aux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);
                        vbef2aux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);
                        vbef3aux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);
                        xBestaux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR) = inicializaSwarm(N_Init/2, N_PAR, X_MAX, X_MIN);  % 初始化新加种群的历史最优位置
                        
                        %%%%%%修改BUG所加内容
                        fitaux{SWARM_ALIVE(1)}=inf*ones(N(SWARM_ALIVE(1)),1);
                        %%%%%%修改BUG所加内容
                        
                        for j=1:N(SWARM_ALIVE(1))   % 计算新加种群的适应度值
                            fitaux{SWARM_ALIVE(1)}(j)=feval(fhd,(xaux{SWARM_ALIVE(1)}(j,:))',varargin{:})-fbias(varargin{:});
                            fitBestaux{SWARM_ALIVE(1)}(j)=fitaux{SWARM_ALIVE(1)}(j);   % 储存新加种群的最优适应度值
                        end
                        [a,b]=min(fitaux{SWARM_ALIVE(1)});
                        gBestaux{SWARM_ALIVE(1)}=xaux{SWARM_ALIVE(1)}(b,:);        % 新加种群的最优位置
                        gbestvalueaux(SWARM_ALIVE(1)) = fitaux{SWARM_ALIVE(1)}(b); % 新加种群的最优位置的适应度值
                        %                                 fprintf('\ncreate swarm %d',i);
                        SC(SWARM_ALIVE(1))=0;
            
                        fit{SWARM_ALIVE(1)}=fitaux{SWARM_ALIVE(1)};
                        x{SWARM_ALIVE(1)}=xaux{SWARM_ALIVE(1)};
                        v{SWARM_ALIVE(1)}=vaux{SWARM_ALIVE(1)};
                        vbef1{SWARM_ALIVE(1)}=vbef1aux{SWARM_ALIVE(1)};
                        vbef2{SWARM_ALIVE(1)}=vbef2aux{SWARM_ALIVE(1)};
                        vbef3{SWARM_ALIVE(1)}=vbef3aux{SWARM_ALIVE(1)};
                        xBest{SWARM_ALIVE(1)}=xBestaux{SWARM_ALIVE(1)};
                        gBest{SWARM_ALIVE(1)}=gBestaux{SWARM_ALIVE(1)};
                        gbestvalue(SWARM_ALIVE(1))=gbestvalueaux(SWARM_ALIVE(1));
                    end
                end
            end
            
        end
    end
    
    fit=fitaux;
    x=xaux;
    v=vaux;
    vbef1=vbef1aux;
    vbef2=vbef2aux;
    vbef3=vbef3aux;
    xBest=xBestaux;
    gBest=gBestaux;
    gbestvalue=gbestvalueaux;
    
    gb(nger)=fitting(nger);
    Convergence_curve(nger)=gb(nger)+fbias(varargin{:});
    nger=nger+1;
end
data = [data;fitting(N_GER)+fbias(varargin{:})];
% toc                      %记录程序运行时间  计时结束
end