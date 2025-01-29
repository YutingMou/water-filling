%TSG 期刊文章 算法2 的仿真程序
clc;
clear;
load('q_nr.mat');
T=168;
%q从17:00算起
q=q_nr';
q(99:168)=q_nr(50:119)';
q(1:98)=cat(2,q_nr(120:168)',q_nr(1:49)');%q从17：00算起

q1=0.2*q;
q2=0.3*q;
q3=0.1*q;
q4=0.25*q;
q5=0.15*q;
Q=[q1;q2;q3;q4;q5];

P=zeros(5,T);%电动汽车负载
p_sum=100;
p1_sum=p_sum/5;
p2_sum=p_sum/5-5;
p3_sum=p_sum/5;
p4_sum=p_sum/5+5;
p5_sum=p_sum/5;
P_sum=[p1_sum,p2_sum,p3_sum,p4_sum,p5_sum];
%P_sum_2=[20 11.592700049728155 14.765953848630190 17.989346136778590 14.193439327493309];
p_max=20;
p1_max=(p_max/5)*ones(1,T);
p1_max(85:168)=zeros(1,84);%假如1号电动车提前走
p2_max=(p_max/5-0.5)*ones(1,T);
p3_max=(p_max/5-0.1)*ones(1,T);
p4_max=(p_max/5+0.6)*ones(1,T);
p5_max=(p_max/5)*ones(1,T);
P_max=[p1_max;p2_max;p3_max;p4_max;p5_max];


for n=1:5
    p_allocated=0;
    a_min=-50;
    a_max=50;
    a=a_min;
while ~(abs(p_allocated-P_sum(n))<0.0000000001)
    a=(a_min+a_max)/2;
for k=1:T %删去P(1,k)则表示此时这个电动车的信息不知道
    
   if n==1
    P(n,k)=a-sum(Q(:,k));
   end
    if n==2
        P(n,k)=a-P(1,k)-sum(Q(:,k));
    end
    
    if n==3 
        P(n,k)=a-sum(Q(:,k))-P(1,k)-P(2,k);
    end
    
    if n==4
        P(n,k)=a-sum(Q(:,k))-P(1,k)-P(2,k)-P(3,k);    
    end
    
    if n==5
        P(n,k)=a-sum(Q(:,k))-P(1,k)-P(2,k)-P(3,k)-P(4,k);    
    end
    
    if P(n,k)<0
        P(n,k)=0;
    else
        if P(n,k)>P_max(n,k)
            P(n,k)=P_max(n,k);
        end
    end
end
p_allocated=sum(P(n,:))/7;%12分钟一个采样周期，因此要除以5
if p_allocated<P_sum(n)
    a_min=a;
else a_max=a;
end
a;
end
end
t=1:T;
figure (1);
plot(t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:),...
t,Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,P(1,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,P(1,:)+P(2,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,P(1,:)+P(2,:)+P(3,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),...
t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:));
legend('P','Q','P+Q','1','2','3','4');
% axis([0,51,-1,30])
figure (2);
plot(t,P(1,:),t,P(2,:),t,P(3,:),t,P(4,:),t,P(5,:));
legend('1','2','3','4','5');
%axis([0,51,-0.5,6])
figure (3);
plot(t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:),':',t,Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),'-.',t,P(1,:)+P(2,:)+P(3,:)+P(4,:)+P(5,:)+Q(1,:)+Q(2,:)+Q(3,:)+Q(4,:)+Q(5,:),'-');
legend('Aggregate PHEV power demand','Aggregate non-PHEV power demand','Total power demand');
%axis([0,51,-1,35])



