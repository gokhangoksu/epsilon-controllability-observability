clc
close all
clear all

%% System Parameters
Lp=-2.93;
Lbeta=-4.75;
Lr=0.78;
gV=0.086;
Ybeta=-0.11;
Nbetadot=0.1;
Np=-0.042;
Nbeta=2.601;
Nr=-0.29;
Lsigmaa=-3.91;
Ysigmar=0.035;
Nsigmar=-2.5335;
Nsigmaa=0.31;

%% System Matrices
A=[0 1 0 0;
   0 Lp Lbeta Lr;
   gV 0 Ybeta -1;
   Nbetadot*gV Np Nbeta+Nbetadot*Ybeta Nr-Nbetadot];
B=[0 0;
   0 Lsigmaa;
   Ysigmar 0;
   Nsigmar+Nbetadot*Ysigmar Nsigmaa];

A={A*0.9;A;A*1.1};
B={B*0.9;B;B*1.1};
K=[-3.4485   -8.0455   -0.2776    1.0903;
    0.0777   -0.6476   -2.7690    4.8456];

%% Time Interval
dt=0.0001;
t1=5;
t=0:dt:t1;

%% System Numerical Solution
x0=[pi/4;0;pi/4;0];
x=zeros(4,t1/dt+1);
u=zeros(2,t1/dt);
x(:,1)=x0;
tic
for k=1:t1/dt
    u(:,k)=K*x(:,k);
    x(:,k+1)=x(:,k)+dt*(A{sigma(t(k))}*x(:,k)+B{sigma(t(k))}*u(:,k));
end
toc

%% Stabilizing State Feedback and State Space Response Plot
plot(t,x,'LineWidth',2)
axis([0 t1 -1 1])
xlabel({'$t$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
legend('x_1(t)','x_2(t)','x_3(t)','x_4(t)','Orientation','horizontal');

%% Final State
x(:,length(t))

%% Sampling epsilon=0.01, n=4
eta=1250; xi=exp(1)/3400;
n=4;

tSampl=[];
P=zeros(round(t1/(3*xi+3/eta)),18);
dt=0.0001;
for m=0:1:round(t1/(3*xi+3/eta))-1
    Pm=[];
    for i=0:1:n-1
        for j=0:1:n-1
            Pm=[Pm (3*m+i)*xi+(3*m+j)/eta];
            tSampl=[tSampl (3*m+i)*xi+(3*m+j)/eta];
        end
    end
    Pm=[Pm (3*m+4)*xi+(3*m)/eta (3*m)*xi+(3*m+4)/eta];
    tSampl=[tSampl (3*m+4)*xi+(3*m)/eta (3*m)*xi+(3*m+4)/eta];
    Pm=sort(Pm);
    P(m+1,:)=Pm;
end
tSampl=round(tSampl,4);
tSampl=sort(tSampl);
tSampl=unique(tSampl);

% Subinterval Check
P415=P(416,:)
P436=P(437,:)
P832=P(833,:)
P874=P(875,:)

%% System Numerical Solution (with Proposed Sampling Schedule)
x0=[pi/4;0;pi/4;0];
x=zeros(4,length(tSampl));
u=zeros(2,length(tSampl));
x(:,1)=x0;
tic
for j=1:length(tSampl)-1
    tau=tSampl(j):dt:tSampl(j+1);
    Intgrnd11=zeros(1,length(tau));
    Intgrnd21=zeros(1,length(tau));
    Intgrnd31=zeros(1,length(tau));
    Intgrnd41=zeros(1,length(tau));
    Intgrnd12=zeros(1,length(tau));
    Intgrnd22=zeros(1,length(tau));
    Intgrnd32=zeros(1,length(tau));
    Intgrnd42=zeros(1,length(tau));
    for k=1:length(tau)
        Int=expm(A{sigma(tau(k))}*(tSampl(j+1)-tau(k)))*B{sigma(tau(k))};
        Intgrnd11(k)=Int(1,1);
        Intgrnd21(k)=Int(2,1);
        Intgrnd31(k)=Int(3,1);
        Intgrnd41(k)=Int(4,1);
        Intgrnd12(k)=Int(1,2);
        Intgrnd22(k)=Int(2,2);
        Intgrnd32(k)=Int(3,2);
        Intgrnd42(k)=Int(4,2);
    end
    Intgrl=[trapz(tau,Intgrnd11) trapz(tau,Intgrnd12);
            trapz(tau,Intgrnd21) trapz(tau,Intgrnd22);
            trapz(tau,Intgrnd31) trapz(tau,Intgrnd32);
            trapz(tau,Intgrnd41) trapz(tau,Intgrnd42)];
    u(:,j)=K*x(:,j);
    x(:,j+1)=expm(A{sigma(tSampl(j))}*(tSampl(j+1)-tSampl(j)))*x(:,j)+Intgrl*u(:,j);
end
toc

%% Stabilizing State Feedback and State Space Response Plot (with Proposed Sampling Schedule)
figure
plot(tSampl,x,'LineWidth',2)
axis([0 t1 -1 1])
xlabel({'$t_j$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
legend('x_{1j}','x_{2j}','x_{3j}','x_{4j}','Orientation','horizontal');

%% Final State and its Length (with Proposed Sampling Schedule)
x(:,length(tSampl))
length(tSampl)