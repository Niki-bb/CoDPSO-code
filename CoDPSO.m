%%%%%%%%%%%%%%%%%2024.2.12 wxb�޸�ע�� ���䣺649713453@qq.com%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;              %������б���
% tic                     %��¼��������ʱ��  ��ʱ��ʼ
% clc;                    %����
function [Convergence_curve,data]=CoDPSO(fhd,pop_size1,dim,Max_iteration,lb,ub,aer,bei,varargin)
fbias=[300, 400, 600, 800, 900, 1800,...
       2000, 2200, 2300, 2400, 2600, 2700];
population=pop_size1;         %Ⱥ�����Ӹ���
N_PAR=dim;               %����ά��
N_GER=Max_iteration;               %����������
nswarms=5;              %Ⱥ����
W=1;
PHI1=1.5;               %ѧϰ����1
PHI2=1.5;               %ѧϰ����2
% [Us]=textread('25DST��ѹ.txt','%.2f');    
% Ut=(Us(2430:12230))';
aerfa=aer;
beita=bei;
Convergence_curve=zeros(1,N_GER); 
data = [];

N_MIN=10;            % 10 ������round(population/2)
N_Init=population;                    % 30
N_MAX=50;                             % ������population*2
MIN_SWARMS=round(nswarms/2);          % round������ʾ��������ȡ�� 3
N_SWARMS=nswarms;                     % 5
MAX_SWARMS=nswarms*2;                 % 10

SWARMS=zeros(MAX_SWARMS,1);           % (10,1)
for i=1:MAX_SWARMS
    if (i<=N_SWARMS)
        SWARMS(i)=1;                  % SWARMS=[1;1;1;1;1;0;0;0;0;0]
    end
end

N=N_Init*ones(MAX_SWARMS,1);          % ά��Ϊ[10,1],����Ԫ��ֵ��ΪN_Init-20
% SCmax = round(N_GER/10);               % ����SCmax��ȡֵ��ο�һЩ���ף��Լ�����һ�� �˴�SCmaxΪ1000
SCmax = 30; 
SC=zeros(MAX_SWARMS,1);               % (10,1)

X_MAX=ub;
X_MIN=lb;
vmax = (X_MAX-X_MIN)/4;
vmin = -vmax;

gbestvalue = 1000000*ones(MAX_SWARMS,1); 
gbestvalueaux=gbestvalue;

fit = cell(MAX_SWARMS,1);       % ��Ϊ(10,1) cell����һ�㱻����Ԫ�����飬����ÿ����Ԫ���Դ��治ͬ���������ͣ���������ֵ���ַ�������Ԫ������ȣ�������ѧ����c������Ľṹ�塣��������Ԫ�������ÿ����Ԫ������Ǿ���
x=cell(MAX_SWARMS,1);
v=cell(MAX_SWARMS,1);
vbef1=cell(MAX_SWARMS,1);
vbef2=cell(MAX_SWARMS,1);
vbef3=cell(MAX_SWARMS,1);
for i=1:N_SWARMS                % i=1:5
    x{i}=zeros(N(i),N_PAR);     % fitness do vector best (20,N_PAR)ȫ����   
    fit{i}(1:N(i),1)=0;         % fitness do vector best 
    v{i}=zeros(N(i),N_PAR);     % inicializa v a NxN_PAR zeros
    vbef1{i}=zeros(N(i),N_PAR);  % ��������Ӧ�ó�ʼ��Ϊ��ͬ�ķ�Χ
    vbef2{i}=zeros(N(i),N_PAR);
    vbef3{i}=zeros(N(i),N_PAR);
end
xaux=x;
fitBest=fit;
vaux=v;
vbef1aux=vbef1;       % ���������ٶ�ֵ�ѽ��з������ٶȸ���
vbef2aux=vbef2;
vbef3aux=vbef3;

xBest=cell(MAX_SWARMS,1);
gBest=cell(MAX_SWARMS,1);
Nkill=zeros(MAX_SWARMS,1);   %(10,1)ȫ����
nger=1;
%%%%%%%%%%%%%%%%%%��ʼ����������λ�ú�����ֵ%%%%%%%%%%%%%%%%%%%
for i=1:N_SWARMS        % i=1:5
    xaux{i}(1:N(i),1:N_PAR)=inicializaSwarm(N(i), N_PAR, X_MAX, X_MIN); %x ser? uma matriz com todos os valores de k  / inicializaSwarm��Ϊ��ʼ������Ⱥ����̫����һ��������ɶ��
    x=xaux;
    for j=1:N(i)
        fit{i}(j)=feval(fhd,(x{i}(j,:))',varargin{:})-fbias(varargin{:});
        fitBest{i}(j)=fit{i}(j);      %INICIALIZAR FITBEST  fitBest����ÿ��Ⱥ����ÿ�����ӵ���Ӧ��ֵ
    end
    
    [a,b]=min(fit{i});
    gBest{i}=x{i}(b,:);           % gBest����ÿ��Ⱥ���е���������ֵ
    gbestvalue(i) = fit{i}(b);    % gbestvalue����ÿ��Ⱥ�����������ӵ���Ӧ��ֵ
    gbestvalueaux(i)=gbestvalue(i);  %
    % Aqui inicializa-se o xBest (Matriz do melhor local para cada part�cula)
    xBest{i}(1:N(i),1:N_PAR) = inicializaSwarm(N(i), N_PAR, X_MAX, X_MIN); %melhor x
end

fitaux=fit;
xBestaux=xBest;
gBestaux=gBest;
fitBestaux=fitBest;  
gb=ones(N_GER,1);                % ��¼����������ÿһ����������Ӧ��ֵ
while (nger<=N_GER)
    for i=1:MAX_SWARMS      % i=1:10
        if (SWARMS(i)==1)   %swarm is alive
            
            xBestaux{i}=xBest{i};
            vbef1aux{i}=vbef1{i};
            vbef2aux{i}=vbef2{i};
            vbef3aux{i}=vbef3{i};
            randnum1 = rand ([N(i), N_PAR]);
            randnum2 = rand ([N(i), N_PAR]);
            
            % �ٶȸ���
            gaux=ones(N(i),1);
            vaux{i} = W*(aerfa*v{i}(1:N(i),:) - (1/2)*(aerfa^2-aerfa-beita^2)*vbef1{i}(1:N(i),:)...
                      + (1/6)*(aerfa^3-3*aerfa^2+aerfa*(2-3*beita^2)+3*beita^2)*vbef2{i}(1:N(i),:)...
                      + (1/24)*(aerfa^4+4*beita^4-4*aerfa^3-6*aerfa^2*beita^2+11*(aerfa^2-beita^2)+18*aerfa*beita^2-6*aerfa)*vbef3{i}(1:N(i),:))...
                      + randnum1.*(PHI1.*(xBest{i}(1:N(i),:)-x{i}(1:N(i),:)))...
                      + randnum2.*(PHI2.*(gaux*gBest{i}-x{i}(1:N(i),:)));
            % �ٶ�����
            for j=1:N_PAR
                for q=1:N(i)
                    vaux{i}(q,j) = ( (vaux{i}(q,j) <= vmin).*vmin ) + ( (vaux{i}(q,j) > vmin).*vaux{i}(q,j) );
                    vaux{i}(q,j) = ( (vaux{i}(q,j) >= vmax).*vmax ) + ( (vaux{i}(q,j) < vmax).*vaux{i}(q,j) );
                end
            end
            vbef3aux{i}=vbef2{i};
            vbef2aux{i}=vbef1{i};
            vbef1aux{i}=v{i};
            % λ�ø���
            xaux{i} = x{i}(1:N(i),:)+vaux{i}(1:N(i),:);

            % limit the position between the mimimum and maximum thresholds
            % λ������
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
                if fitaux{i}(j) < fitBestaux{i}(j)     % fitBestaux�Ǵ洢��ÿ��Ⱥ����ÿ�����ӵ���Ӧ��ֵ���൱�ڳ�ʼ��fitaux��ֵ��fitaux�ǵ�ǰ�ġ�
                   fitBestaux{i}(j) = fitaux{i}(j);
                   xBestaux{i}(j,:) = xaux{i}(j,:);
                end
            end
            
            [a,b]=min(fitaux{i});                 % ��������Ӧ����min(fitBestaux)?
            if (a < gbestvalueaux(i))
                gBestaux{i}=xaux{i}(b,:);         % gBestaux����ÿ��Ⱥ���е���������
                gbestvalueaux(i) = fitaux{i}(b);  % gbestvalueaux����ÿ��Ⱥ���е��������ӵ���Ӧ��ֵ
            end
            [fitting(nger),indexM]=min(gbestvalueaux);  % fitting��������Ⱥ���е�������Ӧ��ֵ
%             gBestT(nger,1:N_PAR) = gBestaux{indexM};    % gBestT��������Ⱥ���о���������Ӧ��ֵ������
            
            if (nger>1)
                if gbestvalueaux(i)>=gbestvalue(i)  %check performance of the swarm  ����Ӧ���ǿ�����Сֵ����ʱ�����ֵ��������Ӧ����С�ڵ��ڣ����������������Ӧ����
%                     if gbestvalueaux(i)>2870      %
%                     DE���ԣ�2870Ϊ�˻��㣬�����㷨�����ж��Ƿ���Ҫ��Ӵ˲���
                        SC(i)=SC(i)+1;
                        if (SC(i)==SCmax)        %reached search limit  SCmax
                            if N(i)>N_MIN        %delete particle  N_MIN
                               [a,b]=max(fitaux{i});    % ��Ϊÿ�γͷ�ֵɾ��һ�����ӣ���������ֻ����һ��b����Ϊ���ܴ��ڶ�����Ӷ����������Ӧ��ֵ��ֻ��Ŀ����ά�Ƚϴ�ʱ������ܱȽ�С��

                                if b==1      % ��������ֵ�ڵ�һ��ʱ����������Сֵ���⣩
                                    xaux{i}=xaux{i}(2:N(i),:);
                                    xBestaux{i}=xBestaux{i}(2:N(i),:);
                                    vaux{i}=vaux{i}(2:N(i),:);
                                    vbef1aux{i}=vbef1aux{i}(2:N(i),:);
                                    vbef2aux{i}=vbef2aux{i}(2:N(i),:);
                                    vbef3aux{i}=vbef3aux{i}(2:N(i),:);
                                end
                                if b==N(i)   % ��������ֵ�����һ��ʱ
                                    xaux{i}=xaux{i}(1:b-1,:);
                                    xBestaux{i}=xBestaux{i}(1:b-1,:);
                                    vaux{i}=vaux{i}(1:b-1,:);
                                    vbef1aux{i}=vbef1aux{i}(1:b-1,:);
                                    vbef2aux{i}=vbef2aux{i}(1:b-1,:);
                                    vbef3aux{i}=vbef3aux{i}(1:b-1,:);
                                end
                                if (b>1 && b<N(i))
                                    xaux{i}=vertcat(xaux{i}(1:(b-1),:),xaux{i}((b+1):N(i),:));    % vercat��ʾ��ֱ��������
                                    xBestaux{i}=vertcat(xBestaux{i}(1:(b-1),:),xBestaux{i}((b+1):N(i),:));
                                    vaux{i}=vertcat(vaux{i}(1:(b-1),:),vaux{i}((b+1):N(i),:));
                                    vbef1aux{i}=vertcat(vbef1aux{i}(1:(b-1),:),vbef1aux{i}((b+1):N(i),:));
                                    vbef2aux{i}=vertcat(vbef2aux{i}(1:(b-1),:),vbef2aux{i}((b+1):N(i),:));
                                    vbef3aux{i}=vertcat(vbef3aux{i}(1:(b-1),:),vbef3aux{i}((b+1):N(i),:));
                                end
                                N(i)=N(i)-1;
                                Nkill(i)=Nkill(i)+1;
                                SC(i)=fix(SCmax*(1-1/(Nkill(i)+1)));     % fix������ʾ���㷽��ȡ��
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
                    if (Nkill(i)>0)           % Nkill��ʼֵΪ0
                        Nkill(i)=Nkill(i)-1;
                    end
                    if (N(i)<N_MAX)     %create particle  N_MAX=40
                        
                        N(i)=N(i)+1;
                        %N_ALIVE=size(xaux{i},1); 
                        xaux{i}(N(i),1:N_PAR)=inicializaSwarm(1, N_PAR, X_MAX, X_MIN);     % Ϊ�¼����ӽ��и�ֵ
                        xBestaux{i}(N(i),1:N_PAR)=inicializaSwarm(1, N_PAR, X_MAX, X_MIN); % ��ʼ���¼����ӵ�����λ��
                        vaux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);     % ��ʼ���¼����ӵ��ٶ�ֵ
                        vbef1aux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);
                        vbef2aux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);
                        vbef3aux{i}(N(i),1:N_PAR)=zeros(1,N_PAR);
                        fitaux{i}(N(i))=feval(fhd,(xaux{i}(N(i),:))',varargin{:})-fbias(varargin{:});
                        fitBestaux{i}(N(i))=fitaux{i}(N(i)); % �����¼����ӵ���Ӧ��ֵ
                        [a,b]=min(fitaux{i});     % �ҳ��¼����Ӻ����Ⱥ�е���С��Ӧ��ֵ����ʱ�¼����ӵ���Ӧ��ֵ�Ѿ�������fitaux��
                        gBestaux{i}=xaux{i}(b,:); % �����¼����Ӻ����Ⱥ������λ��
                        gbestvalueaux(i) = fitaux{i}(b); % �����¼����Ӻ����Ⱥ��������Ӧ��ֵ
                        %                                 fprintf('\ncreate particle in swarm %d',i);
                    end
                end
                if (Nkill(i)==0)&&(N_SWARMS<MAX_SWARMS)     %add new swarm N_SWARMS��ʼֵΪ5��MAX��ʼֵΪ10
                    
                    prob=rand/N_SWARMS;
                    create_swarm=rand;
                    if (create_swarm<prob)
                        SC(i)=0;                        % SC��ʼֵΪ0
                        SWARM_ALIVE=find(SWARMS==0);    % search the swarms that are dead ���б�ʾ����SWARMS�е���0��Ԫ�ص�λ��
                        SWARMS(SWARM_ALIVE(1))=1;       % ����һ������0��SWARMS��ֵΪ1��ԭ��Ϊ0��
                        N_SWARMS=N_SWARMS+1;
                        N(SWARM_ALIVE(1))=N_Init;       % N_InitΪ��Ⱥ��С������20
                        
                        %%% �޸Ĳ���
                        [~,b]=sort(fitaux{i});
                        xaux{SWARM_ALIVE(1)}(1:N_Init/2,:)=xaux{i}(b(1:N_Init/2),:);  % ��ʼ���¼���Ⱥ��λ��
                        vaux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vaux{i}(b(1:N_Init/2),:);       % ��ʼ���¼���Ⱥ���ٶ�ֵ
                        vbef1aux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vbef1aux{i}(b(1:N_Init/2),:);
                        vbef2aux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vbef2aux{i}(b(1:N_Init/2),:);
                        vbef3aux{SWARM_ALIVE(1)}(1:N_Init/2,:)=vbef3aux{i}(b(1:N_Init/2),:);
                        xBestaux{SWARM_ALIVE(1)}(1:N_Init/2,:) = xBestaux{i}(b(1:N_Init/2),:);  % ��ʼ���¼���Ⱥ����ʷ����λ��
                        %%% �޸Ĳ���
                        
                        xaux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=inicializaSwarm(N_Init/2, N_PAR, X_MAX, X_MIN);  % ��ʼ���¼���Ⱥ��λ��
                        vaux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);       % ��ʼ���¼���Ⱥ���ٶ�ֵ
                        vbef1aux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);
                        vbef2aux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);
                        vbef3aux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR)=zeros(N_Init/2,N_PAR);
                        xBestaux{SWARM_ALIVE(1)}(N_Init/2+1:N(SWARM_ALIVE(1)),1:N_PAR) = inicializaSwarm(N_Init/2, N_PAR, X_MAX, X_MIN);  % ��ʼ���¼���Ⱥ����ʷ����λ��
                        
                        %%%%%%�޸�BUG��������
                        fitaux{SWARM_ALIVE(1)}=inf*ones(N(SWARM_ALIVE(1)),1);
                        %%%%%%�޸�BUG��������
                        
                        for j=1:N(SWARM_ALIVE(1))   % �����¼���Ⱥ����Ӧ��ֵ
                            fitaux{SWARM_ALIVE(1)}(j)=feval(fhd,(xaux{SWARM_ALIVE(1)}(j,:))',varargin{:})-fbias(varargin{:});
                            fitBestaux{SWARM_ALIVE(1)}(j)=fitaux{SWARM_ALIVE(1)}(j);   % �����¼���Ⱥ��������Ӧ��ֵ
                        end
                        [a,b]=min(fitaux{SWARM_ALIVE(1)});
                        gBestaux{SWARM_ALIVE(1)}=xaux{SWARM_ALIVE(1)}(b,:);        % �¼���Ⱥ������λ��
                        gbestvalueaux(SWARM_ALIVE(1)) = fitaux{SWARM_ALIVE(1)}(b); % �¼���Ⱥ������λ�õ���Ӧ��ֵ
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
% toc                      %��¼��������ʱ��  ��ʱ����
end