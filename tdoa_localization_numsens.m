%% 2D TDOA Localization example
%% Effects of Number of Sensors
%% Sensor positions are assumed to be perfectly known
clear; close; clc;

xt=[19 24]';   % Target Location (true)
x0=[0 0]';     % Reference Sensor Location (simplifies the equations)
x1=[0 50]';   
x2=[100 0]';
x3=[100 50]';

Ns=15;  % number of sensors

N=0:Ns-4; % sensors other than 4 reference sensors
for nn=1:length(N)
  xsn=100*rand(1,N(nn));
  ysn=50*rand(1,N(nn));
  ps=[xsn;ysn];
  %% Matrix Form of Sensor Locations
  xs=[x0 x1 x2 x3 ps];
  Ns=max(size(xs));
  %%
  c=3e8; % Speed of light (m/s)
  dr = xt.*ones(2,Ns)-xs;
  rs = diag(sqrt(dr'*dr)); % target to sensor range vector
  %%
  MC=500; % Number of Monte Carlo runs
  std=1e-9; % Timing error standard deviation(sec)
  n=std*randn(Ns-1,MC);
  %% Measurement TDOA Generation
  dt=(1/c)*(rs(2:Ns)-rs(1)*ones(Ns-1,1))*ones(1,MC)+n; % TDOAs for each MC RUN

  for k=1:MC
      b=-diag(xs(:,2:end)'*xs(:,2:end))+c^2*dt(:,k).^2;
     %% Linear Solution
      A=-2*[xs(:,2:end);
            c*dt(:,k)']';
      xest1(:,k)=(A'*A)^-1*A'*b;
      er1(:,k)=xt-xest1(1:2,k);
      %% Nonlinear Solution
      d=c*[dt(:,k)];                
      xc=[50 50]'; % initial estimate (center of the area)
      for mm=1:10  % For loop is used instead of checking |dx|<epsilon in while loop
          dy = b + 2*d*sqrt(xc'*xc)+2*xs(:,2:end)'*xc;
          H=-2*[d*(1/sqrt(xc'*xc))*xc(1)+2*xs(1,2:end)' d*(1/sqrt(xc'*xc))*xc(2)+2*xs(2,2:end)'];
          dx=(H'*H)^-1*H'*dy;
          xc=xc+dx;
      end
      xest2(:,k)=xc;
      er2(:,k)=xt-xest2(1:2,k);
  end
  rmse1(nn)=sqrt((1/MC)*sum(sum(er1.^2,2),1));
  rmse2(nn)=sqrt((1/MC)*sum(sum(er2.^2,2),1));
end

figure
semilogy(N,rmse1,'LineWidth',1.5); hold on;
semilogy(N,rmse2,'LineWidth',1.5); hold on; grid on;
legend('LS','NLS');
xlabel('Number of Sensors');
ylabel('2D - RMSE (m)');
title('Estimation Performance vs Number of Sensors');


figure
plot(xt(1),xt(2),'x','linewidth',4,'markersize',20,'color','r'); hold on;
plot(xs(1,:),xs(2,:),'s','markersize',10,'color','b'); hold on; grid on;
plot(xest1(1,:),xest1(2,:),'.','markersize',10,'color','g'); hold on;
plot(xest2(1,:),xest2(2,:),'.','markersize',10,'color','k');  hold on; grid on;
xlabel('X position (m)');
ylabel('Y position (m)');
legend('Target','Sensor','LS','NLS');
str = sprintf('Estimation Results: \\sigma_t=%d ns - %d Sensors',std*1e9,Ns);
title(str);

axes('position',[.65 .175 .25 .25]);
box on % put box around new pair of axes
index= (xest1(1,:) < xt(1)+4*rmse1(end)) & (xest1(1,:) > xt(1)-4*rmse1(end)); % range of estimate indexes around true value
plot(xt(1),xt(2),'x','linewidth',3,'markersize',35,'color','r'); hold on;
plot(xest1(1,index),xest1(2,index),'.','color','g'); hold on;
plot(xest2(1,index),xest2(2,index),'.','color','k'); grid on;% plot on new axes
##axis tight
