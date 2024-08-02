clear
%close
clc
% --- 1. Input Data --------------------------------------
% --- 1.a. Physical --------------------------------------
Re = 40;
% --- 1.b. Numerical -------------------------------------
N = 20; % Cutoff Fourrier mode
c = 0.1; % CFL condition
dt = c * Re/(N^2);
maxIter=1e5; % Maximum number of iterations for the solver
maxRes=1e-3; % Solver residual |Ax-b|<maxRes
sw = 0; % Turbulence model ON/OFF
% --- 2. Previous Calculation and Vector Allocation ------
kvect = (-N:N)'; % allocating vector for Fourrier modes
ki = @(i) i + N+1; % a simple function that links the k values to their corresponding
index
uhat = zeros(length(kvect),maxIter); % allocating uhat for all ks and all timesteps
uhat(:,1) = 1./kvect; % Setting initial values for uhat
uhat(ki(0),1)=0;
conv = zeros(1,N); % Vector for convective term;
dif = zeros(1,N); % Variable for diffusive term
nu_t_star = zeros(1,N);
nu_t = zeros(1,N);
Ek = zeros(N,maxIter); % Allocating energy (only for positive mode)
Ek(1:N,1) = uhat(ki(1:N),1).*conj(uhat(ki(1:N),1)); % Energy for first time-step
% --- 3. Solve using GS ---------------------------------
res = 1 + maxRes;
iter = 1;
t = 0;
% ---- 3.a. The time loop -------------------------------
while res>maxRes && iter<maxIter
uhat(:,iter+1) = uhat(:,iter);
% ---- 3.b. The k loop ------------------------------
for k = 2:N
% ---- 3.b.1 The convective term --------------------conv(k) = 0;
for p =-(N-k):N
q = k-p;
conv(k) = conv(k) + uhat(ki(p),iter+1)*1i*q*uhat(ki(q),iter+1);
end
% ---- 3.b.2 LES SGS model ---------------------------
% ---- 3.b.2 nu_t ------------------------------------
m=2;
Ck = 1; %0.4523;
nu_t_inf = 0.31*(5-m)/(m+1)*sqrt(3-m)*Ck^-(3/2);
nu_t_star(k) = 1+34.5*exp(-3.03*(N/k));
E_end = conj(uhat(ki(N),iter+1))*uhat(ki(N),iter+1);
nu_t(k) = nu_t_inf*(E_end/N)^0.5*nu_t_star(k);
% ---- 3.b.3 Diffusive term ------------------------------------
dif(k) = - (1/Re+sw*nu_t(k))*k^2*uhat(ki(k),iter+1);
% ---- 3.b.4 uhat at next timestep
uhat(ki(k),iter+1) = uhat(ki(k),iter+1) + dt*(-conv(k) + dif(k));
uhat(ki(1),iter+1) = uhat(ki(1),iter);
uhat(ki(-1),iter+1) = uhat(ki(-1),iter);
% ---- 3.b.5 Calculating Ek
Ek(k,iter+1) = uhat(ki(k),iter+1)*conj(uhat(ki(k),iter+1)); % Energy for
first time-step
Ek(1,iter+1) = 1;
end
% ---- 3.b.5 Using the uhat(-k) = conj(uhat(k)) to calculate uhat
% for negative k's
uhat(ki(-N:-1),iter+1) = conj(uhat(ki(N:-1:1),iter+1));
for k = 1:N
% ENERGY EQUATION
%Ek(k,ite+1) = uhat(ki(k),ite+1) * conj(uhat(ki(k),ite+1));
end
% ---- 4. calculating the new time
res = max( abs(uhat(:,iter+1)-uhat(:,iter))/dt )
t = dt + t;
iter = iter + 1;
end
iter = iter - 1;t = 0:dt:t;
%%
figure(1)
loglog(1:N,Ek(:,iter))
%%
x= -10:0.1:10;
% Creating and opening video files for animating P-h and T-s diagrams
v1 = VideoWriter('burgersUnderResolved.avi');
open(v1)
figure(2)
u = zeros(length(x),length(t));
for j = 1:3:iter
for i = 1:length(x)
u(i,j) = sum(uhat(:,j).*exp(1i*kvect*x(i)));
end
subplot(1,2,1)
plot(x,u(:,j))
xlabel('x')
ylabel('u')
title('t = ',num2str(t(j)))
subplot(1,2,2)
loglog(1:N,Ek(:,j))
xlabel('k')
ylabel('E_k')
frame=getframe(gcf);
writeVideo(v1,frame);
end
close(v1);
% closing video files
