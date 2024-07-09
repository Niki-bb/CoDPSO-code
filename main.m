% clear all
% mex cec22_test_func.cpp -DWINDOWS
% 2023.12.19 wxb�޸�
clear all
clc
tic
% core_number=5;            %��Ҫ���õĴ���������
% parpool('local',core_number);
func_num=11;                   % ���Ժ������ CEC2022�Ǵ�1��12
D=20;                         % ά��
Xmin=-100;                    % λ�ñ߽�
Xmax=100;
pop_size=30;                  % ��Ⱥ������
iter_max=D*1000;              % �ܵ������� D*1000
repeat_count = 30;            % �ظ����еĴ��� 
calculate_CoDPSO=zeros(repeat_count,iter_max); 
sendmin=zeros(121,1);
result=zeros(121,4);
MEAN=1./zeros(121,1);
STD=1./zeros(121,1);
MIN=1./zeros(121,1);
fbias1=[300, 400, 600, 800, 900, 1800,...
       2000, 2200, 2300, 2400, 2600, 2700];
Dot_Interval=1;
Dot_Count = iter_max/Dot_Interval; 
CoDPSO_Dot_Data=zeros(1,Dot_Count);  
fhd=str2func('cec22_test_func');   % �������������ʹ��fhd�൱���ڵ���'cec22_test_func'����

for send_i=1
    count=1;
    for aer=0.3
        for bei=0.8
            rand('seed',send_i);
            end_data = [];
            for i = 1:repeat_count
                [cg_curve,data]=CoDPSO(fhd,pop_size,D,iter_max,Xmin,Xmax,aer,bei,func_num);
                end_data = [end_data;data];              % �ظ�30�Σ�ÿһ�ν����������Ӧ��ֵ������end_data��end_data��һֱ�ı�
                calculate_CoDPSO(i,:)=cg_curve;                           % �ظ�30�Σ���i�ε���Ӧ��ֵ���� ÿһ�ε���Ӧ��ֵ���߸�ʽΪ��1��200�� ÿһ����ÿһ��Ѱ�ţ�200�ε��������õ���Ӧ������ ��30��200��
            end
            MEAN1=mean(end_data);
            STD1=std(end_data);
            MIN1=min(end_data);
            if MEAN1<MEAN(count,:) 
                MEAN(count,:)=MEAN1;
                STD(count,:)=STD1;
                MIN(count,:)=MIN1;
                sendmin(count,:)=send_i;
            end
            count=count+1;
        end
    end
end
result(:,1)=MEAN;
result(:,2)=STD;
result(:,3)=MIN;
result(:,4)=sendmin;
% delete(gcp('nocreate'));
toc