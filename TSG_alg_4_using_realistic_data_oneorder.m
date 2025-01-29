% 假如limit是10  用于TSG文章0806版的仿真
% two=1: without protection; two=2: with protection
% pml=1:low level penetration; 2: high level penetration



clc;
clear;
T=1440;
pml=2;%1:low level penetration; 2: high level penetration

%普通负载
load('Q.mat');
%滤波器
windowSize =2;
q1=filter(ones(1,windowSize)/windowSize,1,Q1)';
q2=filter(ones(1,windowSize)/windowSize,1,Q2)';
q3=filter(ones(1,windowSize)/windowSize,1,Q3)';
q4=filter(ones(1,windowSize)/windowSize,1,Q4)';
q5=filter(ones(1,windowSize)/windowSize,1,Q5)';
% q1=Q1';
% q2=Q2';
% q3=Q3';
% q4=Q4';
% q5=Q5';
Q_oo=[Q1';Q2';Q3';Q4';Q5']/2;
Qo=[q1;q2;q3;q4;q5]/2;

%能量需求
energy_need_1=8;
energy_need_2=12;
energy_need_3=30;
energy_need_4=12;
energy_need_5=17;
Energy_need=[energy_need_1,energy_need_2,energy_need_3,energy_need_4,energy_need_5]*pml;

%Q从17:00算起
Q_o=zeros(5,T);
Q=zeros(5,T);
%滤波后
Q(:,1:480)=Qo(:,961:1440);
Q(:,481:1440)=Qo(:,1:960);
%真实值
Q_o(:,1:480)=Q_oo(:,961:1440);
Q_o(:,481:1440)=Q_oo(:,1:960);



P=zeros(5,T);%电动汽车负载
Pt=zeros(5,T);%电动汽车负载
Pout=zeros(5,T);%电动汽车负载

P_sum=Energy_need;
Energy_inter=zeros(5,4);%intermediate energy need



time=[1 720;1 840;61 900;361 840;421 1440];

p1=3.84*pml;
p2=6.6*pml;
p3=10*pml;
p4=3.52*pml;
p5=11.52*pml;
p1_max=p1*zeros(1,T);
p1_max(time(1,1):time(1,2))=p1*ones(1,time(1,2)-time(1,1)+1);

p2_max=p2*zeros(1,T);
p2_max(time(2,1):time(2,2))=p2*ones(1,time(2,2)-time(2,1)+1);

p3_max=p3*zeros(1,T);
p3_max(time(3,1):time(3,2))=p3*ones(1,time(3,2)-time(3,1)+1);

p4_max=p4*zeros(1,T);
p4_max(time(4,1):time(4,2))=p4*ones(1,time(4,2)-time(4,1)+1);
p5_max=p5*zeros(1,T);
p5_max(time(5,1):time(5,2))=p5*ones(1,time(5,2)-time(5,1)+1);

P_max=[p1_max;p2_max;p3_max;p4_max;p5_max];


%根据不同时间段不同的计算顺序
order=cell(4,1);
order{1}=[1 2;2 1]; 
order{2}=[1 2 3;2 3 1;3 1 2];
order{3}=[1 2 4 3;2 4 3 1;3 1 2 4;4 3 1 2];
order{4}=[1 2 4 3 5;2 4 3 5 1;3 5 1 2 4;4 3 5 1 2;5 1 2 4 3];

 
%注意与上面的time吻合起来
order_time=[1 60;61 360;361 420;421 1440];


P_average=zeros(5,T);

%分成num个时间段进行计算
for num=1:length(order)   
   
    baseline=Q(1,:)+Q(2,:)+ Q(3,:)+ Q(4,:)+ Q(5,:);
    
    M=length(order{num});%在这个时间段 PHEV的数量
    exceed=0;
  for two=1:2; 
      baseline=Q(1,:)+Q(2,:)+ Q(3,:)+ Q(4,:)+ Q(5,:);
 for nn=1:M
  
    p_allocated=0;
    a_min=-500;
    a_max=500;
    a=a_min;
    
    while ~(abs(p_allocated- (1-exceed)*P_sum(order{num}(nn)))<0.0000001)
        a=(a_min+a_max)/2;   
            for kk=order_time(num,1):T
                Pt(order{num}(nn),kk)=a-baseline(kk);
                if Pt(order{num}(nn),kk)<0
                    Pt(order{num}(nn),kk)=0;
                else
                    if Pt(order{num}(nn),kk)>P_max(order{num}(nn),kk)
                        Pt(order{num}(nn),kk)=P_max(order{num}(nn),kk);
                    end
                end
            end         
        p_allocated=sum(Pt(order{num}(nn),:))/60;%1分钟一个采样周期，因此要除以60
        if p_allocated<(1-exceed)*P_sum(order{num}(nn))
            a_min=a;
        else a_max=a;
        end
    end

        baseline=baseline+Pt(order{num}(nn),:);  
 end
 
     Px=Pt(1,:)+Pt(2,:)+Pt(3,:)+Pt(4,:)+Pt(5,:)+Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:);
    Py=0;
    for x=order_time(num,1):T
        if Px(x)>10
        Py=Py+Px(x)-10;
        end
    end
Eall=sum(Pt(1,order_time(num,1):T)+Pt(2,order_time(num,1):T)+Pt(3,order_time(num,1):T)+Pt(4,order_time(num,1):T)+Pt(5,order_time(num,1):T))
if Eall>0
exceed=Py/Eall
end

%num
%two
%nn
  end



end




 P=Pt;
t=1:T;

t_peak1= round(energy_need_1/p1*60);
t_peak2= round(energy_need_2/p2*60);
t_peak3= round(energy_need_3/p3*60);
t_peak4= round(energy_need_4/p4*60);
t_peak5= round(energy_need_5/p5*60);
p_peak1=p1*ones(1,t_peak1);
p_peak2=p2*ones(1,t_peak2);
p_peak3=p3*ones(1,t_peak3);
p_peak4=p4*ones(1,t_peak4);
p_peak5=p5*ones(1,t_peak5);
P_lv=Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:);
P_lv(1:t_peak1)=P_lv(1:t_peak1)+p_peak1;
P_lv(1:+t_peak2)=P_lv(1:+t_peak2)+p_peak2;
P_lv(61:60+t_peak3)=P_lv(61:60+t_peak3)+p_peak3;
P_lv(361:360+t_peak4)=P_lv(361:360+t_peak4)+p_peak4;
P_lv(421:420+t_peak5)=P_lv(421:420+t_peak5)+p_peak5;
figure(1);
subplot(2,1,1),plot(t,Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:),t,P_lv);
xlabel('Time (min)');ylabel('Power (kW)');
legend('Aggregate non-PHEV power demand','Total power demand without DSM');%axis([0,1500,0,50])
subplot(2,1,2),plot(t,Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:),t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:)+Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:));
xlabel('Time (min)');ylabel('Power (kW)');
legend('Aggregate non-PHEV power demand','Total power demand with DSM');%axis([0,1500,0,30])



figure(2);
subplot(2,2,1),plot(t,Q_o(1,:));axis([0,1441,-0.5,3.5]);xlabel('Time (min)');ylabel('Power (kW)');title('a');
subplot(2,2,2),plot(t,Q(1,:));axis([0,1441,-0.5,3.5]);xlabel('Time (min)');ylabel('Power (kW)');title('b');
subplot(2,2,3),plot(t,P(1,:));axis([0,1441,-0.5,2]);xlabel('Time (min)');ylabel('Power (kW)');title('c');
subplot(2,2,4),plot(t,P(1,:)+Q_o(1,:));axis([0,1441,-0.5,3.5]);xlabel('Time (min)');ylabel('Power (kW)');title('d');


figure (3);
plot(t,Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:),t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:)+Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:));
xlabel('Time (min)');ylabel('Power (kW)');
legend('Aggregate non-PEV power demand','Total power demand with DSM');%axis([0,1500,0,30])

format long;
pcharge1=sum(P(1,:))/60/Energy_need(1);
pcharge2=sum(P(2,:))/60/Energy_need(2);
pcharge3=sum(P(3,:))/60/Energy_need(3);
pcharge4=sum(P(4,:))/60/Energy_need(4);
pcharge5=sum(P(5,:))/60/Energy_need(5);
pcharge=[pcharge1,pcharge2,pcharge3,pcharge4,pcharge5];
pcharget=1:5;
plot(pcharget,pcharge);


