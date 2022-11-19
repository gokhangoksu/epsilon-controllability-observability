clc
close all
clear all

%% System Parameters
L=10;
c1=1;
c2=2;
c3=.5;
R=1;
t1=120;

%% System Matrices
A0={[0 1/L;-1/c1 -R/L];
    [0 1/L;-1/c2 -R/L];
    [0 1/L;-1/c3 -R/L]};

B={[0;1];[0;1];[0;1]};

C={[0 1/L];[0 1/L];[0 1/L]};

A={A0{1}-B{1}*C{1};
   A0{2}-B{2}*C{2};
   A0{3}-B{3}*C{3}};

%% Sampling epsilon=1, n=2
eta=5; xi=pi/16;
n=2;

tSampl=[];
P=zeros(round(t1*eta/2),6);
for m=0:1:t1*eta/2-1
    Pm=[];
    for i=0:1:n
        for j=0:1:n
            if(i+j<=2)
                Pm=[Pm i*xi+(2*m+j)/eta];
                tSampl=[tSampl (i)*xi+(2*m+j)/eta];
            end
        end
    end
    Pm=sort(Pm);
    P(m+1,:)=Pm;
end
tSampl=sort(tSampl);
tSampl=unique(tSampl);

% Subinterval Check
P48=P(49,:)
P98=P(99,:)
P148=P(149,:)
P198=P(199,:)
P248=P(249,:)
P298=P(299,:)

K=length(tSampl);

%% Output Signal Construction and its Plot
x0=[0.2;0.2];
x=zeros(2,K);
y=zeros(1,K);
x(:,1)=x0;
y(:,1)=C{sigmaEx3(tSampl(1))}*x0;
tic
for k=1:K-1
    x(:,k+1)=expm(A{sigmaEx3(tSampl(k))}*(tSampl(k+1)-tSampl(k)))*x(:,k); 
    y(:,k+1)=C{sigmaEx3(tSampl(k+1))}*x(:,k+1);
end
toc

plot(tSampl,x','LineWidth',2)
xlabel({'$t_k$ (sec)'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
legend('x_{1k}','x_{2k}','Orientation','horizontal');

figure
plot(tSampl,y,'LineWidth',2)
xlabel({'$t_k$ (sec)'},'Interpreter','latex') % R2018a and earlier
ylabel('$y_k$','Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Observability Companion Form Construction (with Proposed Sampling Schedule)
P=size(C{sigmaEx3(tSampl(1))},1);
N=size(C{sigmaEx3(tSampl(1))},2);

Atilde=cell(1,K-1);
Ctilde=cell(1,K);
Ctilde{1}=C{sigmaEx3(tSampl(1))};

tic
for j=1:K-1
    Atilde{j}=expm(A{sigmaEx3(tSampl(j))}*(tSampl(j+1)-tSampl(j)));
    Ctilde{j+1}=C{sigmaEx3(tSampl(j+1))};
end

calC=zeros(K*P,K*N);

calOA=zeros(K*N,K*N);
calOI=zeros(K*N,N);

for j=1:K
    calC(((j-1)*P+1):(j*P),((j-1)*N+1):(j*N))=Ctilde{j};
    calOA(((j-1)*N+1):(j*N),((j-1)*N+1):(j*N))=AtildeProd(Atilde,1,j-1);
    calOI(((j-1)*N+1):(j*N),:)=eye(N);
end
calO=calC*calOA*calOI;
toc

%% Exact Initial State
x0

%% Initial State Estimation (with Proposed Sampling Schedule)
x0est=pinv(calO)*(y')