clc
close all
clear all

%% Time Interval
dt=0.0001;
t1=20;
t=0:dt:t1;

%% Switching Input
SwPer=5;
u=rem(floor(t/SwPer),2)-0.5;

%% System Matrices
A={[0 1;-1 0];[0 1;-1.1 0];[0 1;-1.2 0];[0 1;-1.25 0]};
B=[0; 1];

%% System Numerical Solution
x0=[0.01;0];
x=zeros(2,t1/dt);
x(:,1)=x0;
tic
for k=1:t1/dt
    x(:,k+1)=x(:,k)+dt*(A{sigma(t(k))}*x(:,k)+B*u(k));
end
toc

%% Phase Portrait Plot
subplot(1,2,1)
plot(x(1,:),x(2,:),'LineWidth',2)
xlabel({'$x_1(t)$'},'Interpreter','latex') % R2018a and earlier
ylabel({'$x_2(t)$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Final State
x(:,length(t))

%% Sampling epsilon=1, n=2
eta=5; xi=exp(1)/14;
n=2;

tSampl=[];
P=zeros(round(t1*eta/2),6);
for m=0:1:t1*eta/2-1
    Pm=[];
    for i=0:1:n
        for j=0:1:n
            if(i+j<=2)
                Pm=[Pm i*xi+(2*m+j)/eta];
                tSampl=[tSampl i*xi+(2*m+j)/eta];
            end
        end
    end
    Pm=sort(Pm);
    P(m+1,:)=Pm;
end
tSampl=sort(tSampl);
tSampl=unique(tSampl);

% Subinterval Check
P10=P(11,:)
P23=P(24,:)
P35=P(36,:)

%% Switching Input (with Proposed Sampling Schedule)
u=rem(floor(tSampl/SwPer),2)-0.5;

%% System Numerical Solution (with Proposed Sampling Schedule)
x0=[0.01;0];
x=zeros(2,length(tSampl));
x(:,1)=x0;
tic
for j=1:length(tSampl)-1
    tau=tSampl(j):dt:tSampl(j+1);
    Intgrnd=zeros(2,length(tau));
    for k=1:length(tau)
        Intgrnd(:,k)=expm(A{sigma(tau(k))}*(tSampl(j+1)-tau(k)))*B;
    end
    Intgrl=[trapz(tau,Intgrnd(1,:));
            trapz(tau,Intgrnd(2,:))];
    x(:,j+1)=expm(A{sigma(tSampl(j))}*(tSampl(j+1)-tSampl(j)))*x(:,j)+Intgrl*u(j);
end
toc

%% Phase Portrait Plot (with Proposed Sampling Schedule)
subplot(1,2,2)
plot(x(1,:),x(2,:),'LineWidth',2)
xlabel({'$x_{1j}$'},'Interpreter','latex') % R2018a and earlier
ylabel({'$x_{2j}$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Final State and its Length (with Proposed Sampling Schedule)
x(:,length(tSampl))
length(tSampl)