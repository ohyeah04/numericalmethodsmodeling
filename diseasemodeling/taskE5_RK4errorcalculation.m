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
b=0.8; %infection rate
delta=0.2; %recovery rate

% importing adjacency matrix file from task description; 
% csvread function has parameters (filename,R1,C1,[R1, C1, R2, C2] to set
% up boundaries; here we use 0,0,[0,0,9,9] for N=10:
a=csvread(filename,0,0,[0,0,9,9]); 

% setting up matrix p, is used for calculating sum in eq. 2a, 2b (S and I)
p=[];

%% Using nested loops to get model data

for i = 1:N
    for n = 1:10
        p=sum(a(n,:).*I(i,:)); %p is a matrix which is used to calculate sum in S and I;

        %initializing function handling for S, I, R which will be used in Runge-Kutta Methods:
        fI=@(t,I,R) b*(1-I-R)*p-delta*I;
        fR=@(t,I,R) delta*I;
       
        % setting up time counter:
        t(i+1)=t(i)+h;
        
        
        % Fourth order Runge-Kutta method:
        % Calculating coefficients k1 for Runge-Kutta Mehod:
        k1I(n) = fI(t(i),I(i,n),R(i,n));
        k1R(n) = fR(t(i),I(i,n),R(i,n));
        
        % Calculating coefficients k2 for Runge-Kutta Mehod:
        k2I(n) = fI(t(i)+0.5*h,I(i,n)+0.5*h*k1I(n),R(i,n)+0.5*h*k1R(n));
        k2R(n) = fR(t(i)+0.5*h,I(i,n)+0.5*h*k1I(n),R(i,n)+0.5*h*k1R(n));
        
        % Calculating coefficients k3 for Runge-Kutta Mehod:
        k3I(n) = fI((t(i)+0.5*h),(I(i,n)+0.5*h*k2I(n)),(R(i,n)+0.5*h*k2R(n)));
        k3R(n) = fR((t(i)+0.5*h),(I(i,n)+0.5*h*k2I(n)),(R(i,n)+0.5*h*k2R(n)));
        
        % Calculating coefficients k4 for Runge-Kutta Mehod:
        k4I(n) = fI((t(i)+h),(I(i,n)+k3I(n)*h),(R(i,n)+k3R(n)*h));
        k4R(n) = fR((t(i)+h),(I(i,n)+k3I(n)*h),(R(i,n)+k3R(n)*h));
    
        % Summing all results into final equation:
        I(i+1,n)=I(i,n)+h/6*(k1I(n) + 2*k2I(n) + 2*k3I(n) + k4I(n));
        R(i+1,n)=R(i,n)+h/6*(k1R(n) + 2*k2R(n) + 2*k3R(n) + k4R(n));
						
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

%% Error test
%True values calculated here:
Itrue=4*exp(-(4*t)/5).*(exp((4*t)/5)/3 - 3/10) - exp(-t/5).*((4*exp(t/5))/3 - 13/10);
Rtrue=exp(-t/5).*((4*exp(t/5))/3 - 13/10) - exp(-(4*t)/5).*(exp((4*t)/5)/3 - 3/10);
Strue=1-Itrue-Rtrue;

%Calculating relative error: abs(Approxmate error/True Value)*100
%transposing matrices of true values for convenience; rounding till 10^-4:
Itruetranspose=round(transpose(Itrue),4);
Rtruetranspose=round(transpose(Rtrue),4);
Struetranspose=round(transpose(Strue),4);

%finding absoulute error:
Iabsolute=abs(Ip-Itruetranspose);
Rabsolute=abs(Rp-Rtruetranspose);
Sabsolute=abs(Sp-Struetranspose);

%creating matrix where results will be written
Ierror = [];
Rerror = [];
Serror = [];

%calculating error %
for n=1:1001
Ierror = [Ierror abs(Iabsolute(n)/Itruetranspose(n))*100];
Rerror = [Rerror abs(Rabsolute(n)/Rtruetranspose(n))*100];
Serror = [Serror abs(Sabsolute(n)/Struetranspose(n))*100];
end

%plotting the result

figure('name', 'True Value graphs')
hold on
grid on
plot(t,Itrue, 'LineWidth', 2, 'Color', [1 0 0])
plot(t,Rtrue, 'LineWidth', 2, 'Color', [0 1 0])
plot(t,Strue, 'LineWidth', 2, 'Color', [0 0 1])
legend('I','R','S')
hold off

figure('name', 'Error (%) Graphs')
hold on
grid on
plot(t,Ierror, 'LineWidth', 2, 'Color', [1 0 0])
plot(t,Rerror, 'LineWidth', 2, 'Color', [0 1 0])
plot(t,Serror, 'LineWidth', 2, 'Color', [0 0 1])
legend('I error (%)','R error (%)','S error (%)')
hold off

toc; %elapsed time
