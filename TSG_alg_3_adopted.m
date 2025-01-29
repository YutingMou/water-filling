%TSG 期刊文章 算法3 的仿真程序
%setp1：按照先离开的先灌水计算的顺序算出每个时间段的intermediate target
%step2：按照n种不同顺序计算出不能顺序的功率，求平均值


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
P_sum_2=[20  12.678588048483144 16.267128449455235 19.024960508033466 14.626340120648303];
p_max=20;
p1_max=(p_max/5)*ones(1,T);
p1_max(85:168)=zeros(1,84);%假如1号电动车提前走
p2_max=(p_max/5-0.5)*ones(1,T);
p3_max=(p_max/5-0.1)*ones(1,T);
p4_max=(p_max/5+0.6)*ones(1,T);
p5_max=(p_max/5)*ones(1,T);
P_max=[p1_max;p2_max;p3_max;p4_max;p5_max];
order=[1 2 3 4 5;2 3 4 5 1;3 4 5 1 2; 4 5 1 2 3;5 1 2 3 4];

P_average=zeros(5,T);
M=5;%5个PHEV，所以具有5个不同的顺序
for m=1:M
   baseline=q1+q2+q3+q4+q5; 
for n=1:5
    p_allocated=0;
    a_min=-500;
    a_max=500;
    a=a_min;
while ~(abs(p_allocated-P_sum_2(order(m,n)))<0.0000000001)
    a=(a_min+a_max)/2;
for k=1:84 
    
    P(order(m,n),k)=a-baseline(k);
    if P(order(m,n),k)<0
        P(order(m,n),k)=0;
    else
        if P(order(m,n),k)>P_max(order(m,n),k)
            P(order(m,n),k)=P_max(order(m,n),k);
        end
    end
end
p_allocated=sum(P(order(m,n),:))/7;%12分钟一个采样周期，因此要除以5
if p_allocated<P_sum_2(order(m,n))
    a_min=a;
else a_max=a;
end
a;
end
baseline=baseline+P(order(m,n),:);
P_average(order(m,n),:)=P_average(order(m,n),:)+P(order(m,n),:);
end
end
P(1,:)=P_average(1,:)/M;
P(2,:)=P_average(2,:)/M;
P(3,:)=P_average(3,:)/M;
P(4,:)=P_average(4,:)/M;
P(5,:)=P_average(5,:)/M;
for m=1:M
   baseline=q1+q2+q3+q4+q5; 
for n=1:5
    p_allocated=0;
    a_min=-500;
    a_max=500;
    a=a_min;
while ~(abs(p_allocated-P_sum(order(m,n)))<0.0000000001)
    a=(a_min+a_max)/2;
for k=85:168
    
    P(order(m,n),k)=a-baseline(k);
    if P(order(m,n),k)<0
        P(order(m,n),k)=0;
    else
        if P(order(m,n),k)>P_max(order(m,n),k)
            P(order(m,n),k)=P_max(order(m,n),k);
        end
    end
end
p_allocated=sum(P(order(m,n),:))/7;%12分钟一个采样周期，因此要除以5
if p_allocated<P_sum(order(m,n))
    a_min=a;
else a_max=a;
end
a;
end
baseline=baseline+P(order(m,n),:);
P_average(order(m,n),:)=P_average(order(m,n),:)+P(order(m,n),:);
end
end
P(:,85:168)=P_average(:,85:168)/M;


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
%axis([0,51,-1,30])
figure (2);
plot(t,P(1,:),t,P(2,:),t,P(3,:),t,P(4,:),t,P(5,:));
legend('1','2','3','4','5');
%axis([0,51,-0.5,6]);
%figure (3);
%plot(t,P(1,:)+Q(1,:),t,P(2,:)+Q(2,:),t,P(3,:)+Q(3,:),t,P(4,:)+Q(4,:),t,P(5,:)+Q(5,:));
%legend('1','2','3','4','5');