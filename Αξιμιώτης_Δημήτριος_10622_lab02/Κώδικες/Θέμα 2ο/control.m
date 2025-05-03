clear;

phi_inf=0.01;
phi_0=100;
lamda=100;
rho=10;
k1=10;
k2=20;

t1=[0:0.001:20];
r1=-(pi/1000)*t1.^2+(pi/50)*t1;
rd= @(t) interp1(t1, r1, t, 'linear', 'extrap');


[t,x]=ode15s(@(t,x)odefun(t,x,phi_inf,phi_0,lamda,rho,k1,k2,rd),[0:0.001:20],[0;0]);

x1=real(x(:,1));
phi=(phi_0-phi_inf)*exp(-lamda*t)+phi_inf;

z1=(x1-r1')./phi;
a=-k1*log((1+z1)./(1-z1));

figure();
plot(t,x1,LineWidth=1);
hold on;
plot(t,r1,LineWidth=0.7);
title('Παρακολούθηση τροχιάς r(t)');
ylabel('r(t)');
xlabel('t');
legend({'r(t)','rd(t)'});
grid on;

error=x(:,1)-r1';

figure();
plot(t,error,LineWidth=1);
hold on;
plot(t,-phi,LineWidth=1);
hold on;
plot(t,phi,LineWidth=1);
title('Έλεγχος συνθήκης |r(t)-rd(t)|<φ(t)');
legend({'r(t)-rd(t)','-φ(t)','φ(t)'});
grid on;

rho_t=rho*ones(length(t),1);

error2=x(:,2)-a;

figure();
plot(t,error2,LineWidth=1);
hold on;
plot(t,-rho_t,LineWidth=1);
hold on;
plot(t,rho_t,LineWidth=1);
title('Έλεγχος συνθήκης |dr(t)/dt-α(t)|<ρ');
legend({'dr(t)/dt-α(t)','-ρ','ρ'});
grid on;

function dx=odefun(t,x,phi_inf,phi_0,lamda,rho,k1,k2,rd)

    a1=1.315;
    a2=0.725;
    a3=0.225;
    b1=1.175;

    phi=(phi_0-phi_inf)*exp(-lamda*t)+phi_inf;
    z1=(x(1)-rd(t))/phi;

    T1=log((1+z1)/(1-z1));
    a=-k1*T1;
   
    z2=(x(2)-a)/rho;
    T2=log((1+z2)/(1-z2));

    u=-k2*T2;
    dx=[x(2);-a1*x(2)-a2*sin(x(1))+a3*x(2)^2*sin(2*x(1))+b1*u];
    disp(t);
end