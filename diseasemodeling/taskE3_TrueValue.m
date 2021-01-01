tic;
%%Cleaning previous output and command window
clear; %clearing results
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
R(1,1)=0;

% importing adjacency matrix file from task description; 
% csvread function has parameters (filename,R1,C1,[R1, C1, R2, C2] to set
% up boundaries; here we use 0,0,[0,0,9,9] for N=10:
a=csvread(filename,0,0,[0,0,9,9]); 

% setting up matrix p, is used for calculating sum in eq. 2a, 2b (S and I)
p=[];

%% Using nested loops to get model data

for i = 1:1000
    for n = 1:10
        p = sum(a(n,:).*I(i,:)); %p is a matrix which is used to calculate sum in S and I;
        
        syms Ireal(x) Rreal(x);
        eqn1 = diff(Ireal) == 0.8*(1-Ireal-Rreal)*p-0.2*Ireal;
        eqn2 = diff(Rreal) == 0.2*Ireal;
        eqns = [eqn1 eqn2];
        cond = [Ireal(0)==0.1, Rreal(0)==0];
        [Ireal, Rreal] = dsolve(eqns,cond)
    end
end
toc;

%% From this calcuation we can obtain true value for the given simulation models in manual:
%as a result, we obtain following equations by which true values can be described:

%Ireal = 4*exp(-(4*x)/5)*(exp((4*x)/5)/3 - 3/10) - exp(-x/5)*((4*exp(x/5))/3 - 13/10)
%Rreal = exp(-x/5)*((4*exp(x/5))/3 - 13/10) - exp(-(4*x)/5)*(exp((4*x)/5)/3 - 3/10)

