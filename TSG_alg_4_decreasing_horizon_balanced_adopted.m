%TSG期刊文章 算法4 的仿真程序

clc;
clear;
load('q_nr.mat');
D=2;%天数
T=168*D;

%q=cat(2,q_nr',q_nr');%q从00:00算起

%q从17:00算起
q=q_nr';
q(99:168)=q_nr(50:119)';
q(1:98)=cat(2,q_nr(120:168)',q_nr(1:49)');%q从17：00算起
q=cat(2,q,q);

q1=0.2*q;
q2=0.3*q;
q3=0.1*q;
q4=0.25*q;
q5=0.15*q;
Q=[q1;q2;q3;q4;q5];

P=zeros(5,T);%电动汽车负载
Pt=zeros(5,T);%电动汽车负载
Pout=zeros(5,T);%电动汽车负载
p_sum=100;
p1_sum=p_sum/5;
p2_sum=p_sum/5-5;
p3_sum=p_sum/5;
p4_sum=p_sum/5+5;
p5_sum=p_sum/5;
P_sum=[p1_sum,p2_sum,p3_sum,p4_sum,p5_sum];
Energy_inter=zeros(5,4);%84 98 105时刻的 intermediate energy need
p_max=20;

%time=[1 84;1 98;8 105;15 98;22 112];%entry time and departure time
time=[1 84;1 98;8 105;43 98;50 126];
%time=[1 112;1 112;1 112;1 112;1 112];
p1_max=(p_max/5)*zeros(1,T);
p1_max(time(1,1):time(1,2))=(p_max/5)*ones(1,time(1,2)-time(1,1)+1);
%p1_max(50:127)=zeros(1,12);
%p1_max(1:15)=zeros(1,15);%假如1号电动车后来才来充电
p2_max=(p_max/5-0.5)*zeros(1,T);
p2_max(time(2,1):time(2,2))=(p_max/5-0.5)*ones(1,time(2,2)-time(2,1)+1);
%p2_max(44:50)=zeros(1,7);
p3_max=(p_max/5-0.1)*zeros(1,T);
p3_max(time(3,1):time(3,2))=(p_max/5-0.1)*ones(1,time(3,2)-time(3,1)+1);
%p3_max(1:20)=zeros(1,20);%假如3号电动车后来才来充电
p4_max=(p_max/5+0.6)*zeros(1,T);
p4_max(time(4,1):time(4,2))=(p_max/5-0.1)*ones(1,time(4,2)-time(4,1)+1);
p5_max=(p_max/5)*zeros(1,T);
p5_max(time(5,1):time(5,2))=(p_max/5-0.1)*ones(1,time(5,2)-time(5,1)+1);
%p5_max(39:50)=zeros(1,12);
%p5_max(1:15)=zeros(1,15);%假如5号电动车后来才来充电
P_max=[p1_max;p2_max;p3_max;p4_max;p5_max];
%P_max(:,92:T)=zeros(5,77+168);
%order=[1 2 3 4 5;2 3 4 5 1;3 4 5 1 2; 4 5 1 2 3;5 1 2 3 4];

%根据不同时间段不同的计算顺序
% order=cell(4,1);
% order{1}=[1 2]; 
% order{2}=[1 2 3];
% order{3}=[1 2 4 3];
% order{4}=[1 2 4 3 5];
order=cell(4,1);
order{1}=[1 2;2 1]; 
order{2}=[1 2 3;2 3 1;3 1 2];
order{3}=[1 2 3 4;2 3 4 1;4 1 2 3;3 4 1 2];
order{4}=[1 2 3 4 5;2 3 4 5 1;4 5 1 2 3;3 4 5 1 2;5 1 2 3 4];




%注意与上面的time吻合起来
order_time=[1 7;8 42;43 49;50 126];



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
        p_allocated=sum(Pt(order{num}(nn),:))/7;%12分钟一个采样周期，因此要除以5
        if p_allocated<P_sum(order{num}(nn))
            a_min=a;
        else a_max=a;
        end
    end   
         %
        Energy_inter(order{num}(nn),1)=sum(Pt(order{num}(nn),1:84))/7;
        Energy_inter(order{num}(nn),2)=sum(Pt(order{num}(nn),1:98))/7;
        Energy_inter(order{num}(nn),3)=sum(Pt(order{num}(nn),1:105))/7;
        Energy_inter(order{num}(nn),4)=sum(Pt(order{num}(nn),1:126))/7;
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
        
        if order_time(num,1)+N-1<85
            T_end=84;
            nb=1;
        end
        
        if order_time(num,1)+N-1>84&&order_time(num,1)+N-1<99
            T_end=98;
            nb=2;
        end
         if order_time(num,1)+N-1>98&&order_time(num,1)+N-1<106
            T_end=105;
            nb=3;
        end
        if order_time(num,1)+N-1>105&&order_time(num,1)+N-1<127
            T_end=126;
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
            p_allocated=sum(P(order{num}((mm-1)*M+n),1:T_end))/7;%12分钟一个采样周期，因此要除以5
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
% P(1,j+1:N+j)=P_average(1,j+1:N+j)/M;
% P(2,j+1:N+j)=P_average(2,j+1:N+j)/M;
% P(3,j+1:N+j)=P_average(3,j+1:N+j)/M;
% P(4,j+1:N+j)=P_average(4,j+1:N+j)/M;
% P(5,j+1:N+j)=P_average(5,j+1:N+j)/M;
%Pout(:,j+1:N+j)=P_average(:,j+1:N+j);



%Pout(:,:)=P(:,:)/M;
%P(:,:)=P_average(:,:)/M;
%Pout=P_average/M/15;

t=1:T;

figure (1);
% plot(t,Pout(1,:)+Pout(2,:)+Pout(3,:)+Pout(4,:)+Pout(5,:),...
% t,Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
% t,Pout(1,:)+Pout(2,:)+Pout(3,:)+Pout(4,:)+Pout(5,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
% t,Pout(1,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
% t,Pout(1,:)+Pout(2,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
% t,Pout(1,:)+Pout(2,:)+Pout(3,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
% t,Pout(1,:)+Pout(2,:)+Pout(3,:)+Pout(4,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:));
% legend('P','Q','P+Q','1','2','3','4');
plot(t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:),':',t,Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),'-.',t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),'-');
legend('Aggregate PHEV power demand','Aggregate non-PHEV power demand','Total power demand');
axis([0,T/2,-1,30]);

P(:,:)=Pout(:,:);
figure (2);
plot(t,P(1,:),t,P(2,:),t,P(3,:),t,P(4,:),t,P(5,:));
legend('1','2','3','4','5');
axis([0,T/2,-1,6]);

figure (3);
plot(t,Pt(1,:),t,Pt(2,:),t,Pt(3,:),t,Pt(4,:),t,Pt(5,:));
legend('1','2','3','4','5');

figure (4);
plot(t,Pt(1,:)+Pt(2,:)+Pt(3,:)+Pt(4,:)+Pt(5,:),...
t,Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,Pt(1,:)+Pt(2,:)+Pt(3,:)+Pt(4,:)+Pt(5,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,Pt(1,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,Pt(1,:)+Pt(2,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,Pt(1,:)+Pt(2,:)+Pt(3,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,Pt(1,:)+Pt(2,:)+Pt(3,:)+Pt(4,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:));
legend('P','Q','P+Q','1','2','3','4');
axis([0,T/2,-1,30]);