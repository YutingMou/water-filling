%仿真所需时间较长，耐心等待。
clc;
clear;
T=1440;

%普通负载
load('Q.mat');
%滤波器
windowSize = 10;
q1=filter(ones(1,windowSize)/windowSize,1,Q1)';
q2=filter(ones(1,windowSize)/windowSize,1,Q2)';
q3=filter(ones(1,windowSize)/windowSize,1,Q3)';
q4=filter(ones(1,windowSize)/windowSize,1,Q4)';
q5=filter(ones(1,windowSize)/windowSize,1,Q5)';
%q1=Q1';
%q2=Q2';
%q3=Q3';
%q4=Q4';
%q5=Q5';
Q_oo=[Q1';Q2';Q3';Q4';Q5'];
Qo=[q1;q2;q3;q4;q5];

%能量需求
energy_need_1=8;
energy_need_2=12;
energy_need_3=30;
energy_need_4=12;
energy_need_5=17;
Energy_need=[energy_need_1,energy_need_2,energy_need_3,energy_need_4,energy_need_5];

%Q从17:00算起
Q_o=zeros(5,T);
Q=zeros(5,T);
Q_o(:,1:480)=Q_oo(:,961:1440);
Q_o(:,481:1440)=Q_oo(:,1:960);
Q(:,1:480)=Qo(:,961:1440);
Q(:,481:1440)=Qo(:,1:960);


P=zeros(5,T);%电动汽车负载
Pt=zeros(5,T);%电动汽车负载
Pout=zeros(5,T);%电动汽车负载

P_sum=Energy_need;
Energy_inter=zeros(5,4);%intermediate energy need



time=[1 720;1 840;61 900;361 840;421 1440];
p1=3.84;
p2=6.6;
p3=10;
p4=3.52;
p5=11.52;
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
order{3}=[1 2 3 4;2 3 4 1;4 1 2 3;3 4 1 2];
order{4}=[1 2 3 4 5;2 3 4 5 1;4 5 1 2 3;3 4 5 1 2;5 1 2 3 4];


%注意与上面的time吻合起来
order_time=[1 60;61 360;361 420;421 1440];


P_average=zeros(5,T);

%分成num个时间段进行计算
for num=1:length(order)   
   
for N=1:order_time(num,2)-order_time(num,1)+1   
    baseline=q1+q2+q3+q4+q5; 
    
    M=length(order{num});%在这个时间段 PHEV的数量
    %%%  get intermediate target according the First Circular Order
 for nn=1:M
    p_allocated=0;
    a_min=-500;
    a_max=500;
    a=a_min;
    while ~(abs(p_allocated- P_sum(order{num}(nn)))<0.0000001)
        a=(a_min+a_max)/2;   
            for kk=order_time(num,1)+N-1:T
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
        if p_allocated<P_sum(order{num}(nn))
            a_min=a;
        else a_max=a;
        end
    end   
         %
        Energy_inter(order{num}(nn),1)=sum(Pt(order{num}(nn),1:720))/60;
        Energy_inter(order{num}(nn),2)=sum(Pt(order{num}(nn),1:840))/60;
        Energy_inter(order{num}(nn),3)=sum(Pt(order{num}(nn),1:900))/60;
        Energy_inter(order{num}(nn),4)=sum(Pt(order{num}(nn),1:1440))/60;
        baseline=baseline+Pt(order{num}(nn),:);  
 end

 %按照M种不同的计算顺序
 P_average=zeros(5,T);
for mm=1:M    
    baseline=q1+q2+q3+q4+q5; 
    
    for n=1:M
        p_allocated=0;
        a_min=-500;
        a_max=500;
        a=a_min;
        
        if order_time(num,1)+N-1<721
            T_end=720;
            nb=1;
        end
        
        if order_time(num,1)+N-1>720&&order_time(num,1)+N-1<841
            T_end=840;
            nb=2;
        end
         if order_time(num,1)+N-1>840&&order_time(num,1)+N-1<901
            T_end=900;
            nb=3;
        end
        if order_time(num,1)+N-1>900&&order_time(num,1)+N-1<1441
            T_end=1440;
            nb=4;
        end        
        
        
        while ~(abs(p_allocated- Energy_inter(order{num}((mm-1)*M+n),nb))<0.0000001)
            a=(a_min+a_max)/2;   
                for k=order_time(num,1)+N-1:T_end
                    P(order{num}((mm-1)*M+n),k)=a-baseline(k);
                    if P(order{num}((mm-1)*M+n),k)<0
                        P(order{num}((mm-1)*M+n),k)=0;
                    else
                        if P(order{num}((mm-1)*M+n),k)>P_max(order{num}((mm-1)*M+n),k)
                            P(order{num}((mm-1)*M+n),k)=P_max(order{num}((mm-1)*M+n),k);
                        end
                    end
                end         
            p_allocated=sum(P(order{num}((mm-1)*M+n),1:T_end))/60;%1分钟一个采样周期，因此要除以60
            if p_allocated<Energy_inter(order{num}((mm-1)*M+n),nb)
                a_min=a;
            else a_max=a;
            end
        end
               
        
        baseline=baseline+P(order{num}((mm-1)*M+n),:);
        %num
        %mm
        %n
        %order{num}((mm-1)*M+n)
    end
    P_average(:,order_time(num,1)+N-1:T)=P_average(:,order_time(num,1)+N-1:T)+P(:,order_time(num,1)+N-1:T);    
end 
 
 Pout(:,order_time(num,1)+N-1:T)=P_average(:,order_time(num,1)+N-1:T)/M;
 P(:,order_time(num,1)+N-1:T)=P_average(:,order_time(num,1)+N-1:T)/M;
%  num
%  n
%  mm
%  a
%  kk
%  k
%  Energy_inter
%  p_allocated
%  Energy_inter(order{num}((mm-1)*mm+n),1)
end

end

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
figure(3);
subplot(2,1,1),plot(t,Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:),t,P_lv);
xlabel('Time (min)');ylabel('Power (kW)');
legend('Aggregate non-PHEV power demand','Total power demand without DSM');axis([0,1500,0,50])
subplot(2,1,2),plot(t,Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:),t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:)+Q_o(1,:)+Q_o(2,:)+Q_o(3,:)+Q_o(4,:)+Q_o(5,:));
xlabel('Time (min)');ylabel('Power (kW)');
legend('Aggregate non-PHEV power demand','Total power demand with DSM');axis([0,1500,0,30])



figure(7);
subplot(2,2,1),plot(t,Q_o(1,:));axis([0,1441,-0.5,7]);xlabel('Time (min)');ylabel('Power (kW)');title('a');
subplot(2,2,2),plot(t,Q(1,:));axis([0,1441,-0.5,7]);xlabel('Time (min)');ylabel('Power (kW)');title('b');
subplot(2,2,3),plot(t,P(1,:));axis([0,1441,-0.5,1.5]);xlabel('Time (min)');ylabel('Power (kW)');title('c');
subplot(2,2,4),plot(t,P(1,:)+Q_o(1,:));axis([0,1441,-0.5,7]);xlabel('Time (min)');ylabel('Power (kW)');title('d');