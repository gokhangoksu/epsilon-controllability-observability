clc
close all
clear all

%% System Parameters
Mass=400;
mass=50;
kStiff=2*10^4;
ktStiff=2.5*10^5;
cMin=3*10^2;
cMax=3.9*10^2;
t1=1;

%% System Matrices
A={[0 1 0 0;
    -kStiff/Mass -cMin/Mass kStiff/Mass cMin/Mass;
    0 0 0 1;
    kStiff/mass cMin/mass -(kStiff+ktStiff)/mass -cMin/mass];
   [0 1 0 0;
    -kStiff/Mass -cMax/Mass kStiff/Mass cMax/Mass;
    0 0 0 1;
    kStiff/mass cMax/mass -(kStiff+ktStiff)/mass -cMax/mass]};

B={[0;0;0;ktStiff/mass];[0;0;0;ktStiff/mass]};

C={[-kStiff/Mass -cMin/Mass kStiff/Mass cMin/Mass;
    0 1 0 -1];
    [-kStiff/Mass -cMax/Mass kStiff/Mass cMax/Mass;
    0 1 0 -1]};

%% Sampling epsilon=0.01, n=4
eta=1400; xi=pi/4400;
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
P68=P(69,:)
P103=P(104,:)
P162=P(163,:)
P197=P(198,:)

K=length(tSampl);

%% Output Signal and Output Feedback Construction and their Plot
x0=[0;0;0;0];
x=zeros(4,K);
y=zeros(2,K);
u=zeros(1,K);
x(:,1)=x0;
y(:,1)=C{sigmaEx4(tSampl(1))}*x0;
tic
for k=1:K-1
   tau=tSampl(k):dt:tSampl(k+1);
   Intgrnd11=zeros(1,length(tau));
   Intgrnd21=zeros(1,length(tau));
   Intgrnd31=zeros(1,length(tau));
   Intgrnd41=zeros(1,length(tau));
   for kk=1:length(tau)
       Int=expm(A{sigmaEx4(tau(kk))}*(tSampl(kk+1)-tau(kk)))*B{sigmaEx4(tau(kk))};
       Intgrnd11(kk)=Int(1,1);
       Intgrnd21(kk)=Int(2,1);
       Intgrnd31(kk)=Int(3,1);
       Intgrnd41(kk)=Int(4,1);
   end
   u(:,k)=w(tSampl(k));
   x(:,k+1)=expm(A{sigmaEx4(tSampl(k))}*(tSampl(k+1)-tSampl(k)))*x(:,k)...
       +[trapz(tau,Intgrnd11);
         trapz(tau,Intgrnd21);
         trapz(tau,Intgrnd31);
         trapz(tau,Intgrnd41)].*u(:,k); 
   y(:,k+1)=C{sigmaEx4(tSampl(k+1))}*x(:,k+1);
end
toc

plot(tSampl,x','LineWidth',2)
xlabel({'$t_k$ (sec)'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
legend('x_{1k}','x_{2k}','x_{3k}','x_{4k}','Orientation','horizontal');

figure
plot(tSampl,y','LineWidth',2)
xlabel({'$t_k$ (sec)'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
legend('y_{1k}','y_{2k}','Orientation','horizontal');

figure
plot(tSampl,u','LineWidth',2)
axis([0 max(tSampl) 0 0.6])
xlabel({'$t_k$ (sec)'},'Interpreter','latex') % R2018a and earlier
ylabel({'$u_{k}$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Observability Companion Form Construction (with Proposed Sampling Schedule)
P=size(C{sigmaEx4(tSampl(1))},1);
N=size(C{sigmaEx4(tSampl(1))},2);
M=size(B{sigmaEx4(tSampl(1))},2);

Atilde=cell(1,K-1);
Btilde=cell(1,K-1);
Ctilde=cell(1,K);
Ctilde{1}=C{sigmaEx4(tSampl(1))};

tic
for j=1:K-1
   tau=tSampl(j):dt:tSampl(j+1);
   Intgrnd11=zeros(1,length(tau));
   Intgrnd21=zeros(1,length(tau));
   Intgrnd31=zeros(1,length(tau));
   Intgrnd41=zeros(1,length(tau));
   for k=1:length(tau)
       Int=expm(A{sigmaEx4(tau(k))}*(tSampl(k+1)-tau(k)))*B{sigmaEx4(tau(k))};
       Intgrnd11(k)=Int(1,1);
       Intgrnd21(k)=Int(2,1);
       Intgrnd31(k)=Int(3,1);
       Intgrnd41(k)=Int(4,1);
   end
   Atilde{j}=expm(A{sigmaEx4(tSampl(j))}*(tSampl(j+1)-tSampl(j)));
   Btilde{j}=[trapz(tau,Intgrnd11);
              trapz(tau,Intgrnd21);
              trapz(tau,Intgrnd31);
              trapz(tau,Intgrnd41)];
    Ctilde{j+1}=C{sigmaEx4(tSampl(j+1))};
end

calC=zeros(K*P,K*N);

calOA=zeros(K*N,K*N);
calOI=zeros(K*N,N);

calTA=zeros(K*N,K*N);
calTB=zeros(K*N,K*M);

for j=1:K
   calC(((j-1)*P+1):(j*P),((j-1)*N+1):(j*N))=Ctilde{j};
   calOA(((j-1)*N+1):(j*N),((j-1)*N+1):(j*N))=AtildeProd(Atilde,1,j-1);
   calOI(((j-1)*N+1):(j*N),:)=eye(N);
   for jj=2:(j-2)
       calTA(((jj-1)*N+1):(jj*N),((jj-1)*N+1):(jj*N))=AtildeProd(Atilde,jj+1,jj-1);
   end
   if j~=K
       calTB(((j-1)*N+1):(j*N),((j-1)*M+1):(j*M))=Btilde{j};
   end
end
calO=calC*calOA*calOI;
calT=calC*calTA*calTB;
toc

%% Exact Initial State
x0

%% Initial State Estimation (with Proposed Sampling Schedule)
yNew=zeros(1,K*P);
yNew(1:2:K*P)=y(1,:);
yNew(2:2:K*P)=y(2,:);
x0est=pinv(calO)*(yNew'-calT*u')