

% Guided Object and its Moving Guide (pursuer and evador)

clear all
close all
clc

%--------------------Parameters

N = 250;

T = 1;  %7

F = [1 T 0 0;...
    0 1 0 0;...
    0 0 1 T;...
    0 0 0 1];

qx = 0.005;
qy = 0.005;

Q = [qx*T^3/3 qx*T^2/2 0 0;...
    qx*T^2/2 qx*T 0 0;...
    0 0 qy*T^3/3 qy*T^2/2;...
    0 0 qy*T^2/2 qy*T];

D = chol(Q,'lower');

%================================= 

Y1_m_r = [5000;15;2000;10];
Y1_Cov_r = diag([1000;10;1000;10]);

Yr = zeros(4,N);

Dyn = chol(Y1_Cov_r,'lower');
Yr(:,1) = Y1_m_r + Dyn*randn(4,1);

%--------------------Moving Guide

Dyn = chol(Y1_Cov_r);
Yr(:,1) = Y1_m_r + Dyn*randn(4,1);

for k=2:N

     if k>1 && k<100


    v = randn(4,1);
    Yr(:,k) = F*Yr(:,k-1) + D*v;

%     Yr(:,k) = Yr(:,k-1);

    
    elseif k>99 && k<115

            v = randn(4,1);
    Yr(:,k) = F*Yr(:,k-1) + [-2;-1;3;2] + D*v;

    elseif k>114 % && k<250

        v = randn(4,1);
        Yr(:,k) = F*Yr(:,k-1) + D*v;
% 
%     elseif k>59 && k<75
% 
%         v = randn(4,1);
%         Yr(:,k) = F*Yr(:,k-1) + [-15;-5;3;2] + D*v;
%         
%     elseif k>74
%         
%         v = randn(4,1);
%         Yr(:,k) = F*Yr(:,k-1) + D*v;
% 
    end    

end

%--------------------------------

%---------------

Xr = zeros(4,N);

X1_m_r = [5000;10;3000;-10];
X1_Cov_r = diag([1000;10;1000;10]);

%----------

Dxn = chol(X1_Cov_r,'lower');
Xr(:,1) = X1_m_r + Dxn*randn(4,1);


%+++++++++++++++++++++++++++++++++++++ Guided Object

N1 = N;

for k=2:N-1


    %-----------CN|k

    CNk = zeros(4,4);

    for ii=0:N1-2-1

        CNk = CNk + F^ii*Q*(F^ii)';

    end

    %----------

    Covp = Q - Q*(F^(N1-2))'/(CNk + F^(N1-2)*Q*(F^(N1-2))')*F^(N1-2)*Q;

    %----------

    Mp = F*Xr(:,k-1) + Covp*(F^(N1-2))'/(CNk)*(Yr(:,k-1) - F^(N1-2+1)*Xr(:,k-1));

    %----------

    Dcovp = chol(Covp,'lower');

    %----------

    Dc = Dcovp*randn(4,1);
    Xr(:,k) = Mp + Dc;

    %..........

    N1 = N1-1;
    
end

Xr(:,N) = Yr(:,N);

D_EE_p = (sqrt((Xr(1,:)-Yr(1,:)).^2 + (Xr(3,:)-Yr(3,:)).^2));
D_EE_v = (sqrt((Xr(2,:)-Yr(2,:)).^2 + (Xr(4,:)-Yr(4,:)).^2));

%--------------------


figure(1)
hold on
plot(Yr(1,:),Yr(3,:),'.r')
hold on
plot(Xr(1,:),Xr(3,:),'.')


figure(2)
hold on
plot(D_EE_p)

figure(3)
hold on
plot(D_EE_v)


