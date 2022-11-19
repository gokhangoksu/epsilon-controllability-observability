clc
close all
clear all

%% Time Interval
dt=0.0001;
t1=20;
t=0:dt:t1;

%% System Matrices
A={[0 1;-1 0];[0 1;-1.1 0];[0 1;-1.2 0];[0 1;-1.25 0]};
B=[0; 1];
K=[0.5345 -1.1931];     % Feedback Gain

%% System Numerical Solution
x0=[1;0];
x=zeros(2,t1/dt);
u=zeros(1,t1/dt);
x(:,1)=x0;
tic
for k=1:t1/dt
    x(:,k+1)=x(:,k)+dt*((A{sigma(t(k))}+B*K)*x(:,k));
    u(k)=K*x(:,k);
end
toc

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
P10=P(11,:)
P23=P(24,:)
P35=P(36,:)

%% System Numerical Solution (with Proposed Sampling Schedule)
x0=[1;0];
x=zeros(2,length(tSampl));
u=zeros(1,length(tSampl));
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
    x(:,j+1)=expm(A{sigma(tSampl(j))}*(tSampl(j+1)-tSampl(j)))*x(:,j)+Intgrl*K*x(:,j);
    u(j)=K*x(:,j);
end
toc

%% Stabilizing State Feedback and State Space Response Plot (with Proposed Sampling Schedule)
subplot(1,2,1)
plot(tSampl(1:length(u)),u,'LineWidth',2)
axis([0 t1 -2 2])
xlabel({'$t_j$'},'Interpreter','latex') % R2018a and earlier
ylabel({'$u_j$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

subplot(1,2,2)
plot(tSampl,x,'LineWidth',2)
xlabel({'$t_j$'},'Interpreter','latex') % R2018a and earlier
set(gca,'fontsize',24)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
legend('x_{1j}','x_{2j}');

%% Final State and its Length (with Proposed Sampling Schedule)
x(:,length(tSampl))
length(tSampl)