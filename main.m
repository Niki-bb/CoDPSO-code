% clear all
% mex cec22_test_func.cpp -DWINDOWS
% 2023.12.19 wxb修改
clear all
clc
tic
% core_number=5;            %想要调用的处理器个数
% parpool('local',core_number);
func_num=11;                   % 测试函数序号 CEC2022是从1到12
D=20;                         % 维度
Xmin=-100;                    % 位置边界
Xmax=100;
pop_size=30;                  % 种群粒子数
iter_max=D*1000;              % 总迭代次数 D*1000
repeat_count = 30;            % 重复运行的次数 
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
fhd=str2func('cec22_test_func');   % 函数句柄，后续使用fhd相当于在调用'cec22_test_func'函数

for send_i=1
    count=1;
    for aer=0.3
        for bei=0.8
            rand('seed',send_i);
            end_data = [];
            for i = 1:repeat_count
                [cg_curve,data]=CoDPSO(fhd,pop_size,D,iter_max,Xmin,Xmax,aer,bei,func_num);
                end_data = [end_data;data];              % 重复30次，每一次结束后最佳适应度值储存在end_data，end_data在一直改变
                calculate_CoDPSO(i,:)=cg_curve;                           % 重复30次，第i次的适应度值曲线 每一次的适应度值曲线格式为（1，200） 每一行是每一次寻优（200次迭代）所得的适应度曲线 （30，200）
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