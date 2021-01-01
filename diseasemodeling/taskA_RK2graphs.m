tic; %tic-toc for function exectution time calculation;

%%Cleaning previous output and command window
clear; %clearing results
clf; %clearing figure
clc; %clearing command window

% adjacency matrix file:
filename = 'AM_N=50.csv';

%% Initial conditions:

% controlling stepsize and maximum time of the model (by default 100 days)
h=0.1;
tfinal=100;
N=ceil(tfinal/h);

% creating matrices to store data of healthy/infected/recovered people
I = zeros(N,10);
R = zeros(N,10);
%S = ones(N,10);

% initial conditions: note that matlab's initial conditions are set from 1,
% not from zero (eg. t(1)=0);
t(1)=0;
S(1,1)=0;
I(1,2)=1;
I(4,2)=1;
R(1,1)=0;

%population parameters:
b=0.3; %infection rate
delta=0.1; %recovery rate

% importing adjacency matrix file from task description; 
% csvread function has parameters (filename,R1,C1,[R1, C1, R2, C2] to set
% up boundaries; here we use 0,0,[0,0,9,9] for N=10:
a=csvread(filename,0,0,[0,0,9,9]); 

% setting up matrix p, is used for calculating sum in eq. 2a, 2b (S and I)
p=[];

%% Using nested loops to get model data

for i = 1:N
    for n = 1:10
        p = sum(a(n,:).*I(i,:)); %p is a matrix which is used to calculate sum in S and I;

        %initializing function handling for I, R which will be used in Runge-Kutta Methods:
        %overall, we have 3 equations(S, I and R), however, we can solve
        %two and find S, because P(S)+P(I)+P(R)=1; where P - probability;
        fI=@(t,I,R) b*(1-I-R)*p-delta*I;
        fR=@(t,I,R) delta*I;
       
        % setting up time counter:
        t(i+1)=t(i)+h;
        
        % Second order Runge-Kutta method:
        % Calculating coefficients k1 for Runge-Kutta Mehod:
        k1I(n) = fI(t(i),I(i,n),R(i,n));
        k1R(n) = fR(t(i),I(i,n),R(i,n));
        
        % Calculating coefficients k2 for Runge-Kutta Method:
        k2I(n) = fI(t(i)+h/2, I(i,n)+h/2*k1I(n), R(i,n)+h/2*k1R(n));
        k2R(n) = fR(t(i)+h/2, I(i,n)+h/2*k1I(n) ,R(i,n)+h/2*k1R(n));
            
        % Summing all results into final equation:
        I(i+1,n)=I(i,n)+h/2*(k1I(n) + k2I(n));
        R(i+1,n)=R(i,n)+h/2*(k1R(n) + k2R(n));
     end
end

%now, finding S, which is:
S = 1-I-R;

%Now, summing all probabilities to plot probabilites for all i's:
%Sp, Ip, Rp will be total probability, stating them:
Sp = 0;
Ip = 0;
Rp = 0;

%Calculating:
for n=1:10
    Sp=Sp+S(:,n);
    Ip=Ip+I(:,n);
    Rp=Rp+R(:,n);
end

%now we have a sum of 10 probabilities, dividing it by 10 to scale by 1.0:
Sp = Sp/10;
Ip = Ip/10;
Rp = Rp/10;


%% Plotting resultant data
% to set up parameters of graph, using hold command
hold on

%plotting S,I and R of node 10 vs. time, then setting up line width as 2 and
%setting up graphs' colors 
plot(t,Sp, 'LineWidth', 2, 'Color', [0 0 1])
plot(t,Ip, 'LineWidth', 2, 'Color', [1 0 0])
plot(t,Rp, 'LineWidth', 2, 'Color', [0 1 0])

%legend of the graph for better understanding and navigation
legend('S - node being healthy and not infected','I - node being infected by the virus','R - node recovering from the virus')
grid on %turning on the grid
xlim auto  %automatic scaling of x axis limits
title('Modeling result for Second order Runge-Kutta Method') %setting up the title of the plot
xlabel('time (t) days') %x-axis label 
ylabel('Density') %y-axis lablel
%closing hold 
hold off


toc; %elapsed time