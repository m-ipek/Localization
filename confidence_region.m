%% Confidence Region Example for 2 Parameter Curve Fitting Problem
clear; close; clc;
%% y = a*t + b --> Curve Fitting Problem
%% x =[a b]' --> Parameter Vector
t=0:0.1:5; t=t'; % time vector
L=length(t);
xt=[3;2]; % True polynomial coefficients
H=[t ones(L,1)];
yt=H*xt;

std=0.5; % Noise standard deviation
R=std^2*eye(L);

##plot(t,yt); hold on;
##plot(t,y); grid on;

%% ESTIMATOR
MC=100;
for k=1:MC
  n=std*randn(L,1); % zero mean Gaussian noise with sigma=std.
  y=yt+n; % noisy measurement vector
  x(:,k)=(H'*R^-1*H)^-1*H'*R^-1*y; % weighted LS solution
end 
mx=(1/MC)*sum(x,2);
er=x - mx.*ones(2,MC);
%% er=x - xt.*ones(2,MC);
C=(1/MC)*er*er'; % Sampled Estimation Error Covariance Matrix 
Ct=(H'*R^-1*H)^-1; % Theoretical Estimation Error Covariance Matrix

K=5.99; % 1.38 (for p=0.5) AND 5.99 (for p = 0.95): From Chi2 Tables or Chi2 CDF

%% Ellipse Parameters
[V,D] = eig(C); % SVD
d1=D(1,1); % first eigenvalue
d2=D(2,2); % second eigenvalue
v1=V(:,1)/norm(V(:,1));
v2=V(:,2)/norm(V(:,2));
P=[v1 v2]; % Transformation Matrix;

% x^2/ax^2+y^2/ay^2=K --> ELLIPSE EQUATION
ex=-sqrt(d1*K):1e-6:sqrt(d1*K);
ey=d2^0.5*sqrt(K-(ex/d1^0.5).^2);
ney=-ey;

elxy=P*[ex;ey];   % concave part
elxyn=P*[ex;ney]; % convex part of the ellipse

elxy=elxy+mx; % mean shifted to the estimate
elxyn=elxyn+mx;

plot(xt(1),xt(2),'x','markersize',10,'linewidth',3,'color','g'); hold on;
plot(x(1,:),x(2,:),'*','markersize',5,'color','b'); hold on;
plot(elxy(1,:),elxy(2,:),'-','markersize',2,'color','r'); hold on;
plot(elxyn(1,:),elxyn(2,:),'-','markersize',2,'color','r'); hold on; grid on;
legend('True','Estimates','95% Confidence');
str=sprintf('Curve Fitting Performance (%d Monte Carlo) --> a * t + b',MC);
title(str);
xlabel('a'); ylabel('b');



