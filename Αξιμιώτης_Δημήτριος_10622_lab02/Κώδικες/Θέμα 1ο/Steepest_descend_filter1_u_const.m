clear;

am=1;
sys1=tf([-1 0],[1 am]);
sys2=tf([0 -1],[1 am]);
sys3=tf([0 1],[1 am]);
sys4=tf([0 1],[1 0]);

t=0 :0.0001 :20;
u=2.5*ones(length(t),1);
t_data=t;
u_input = @(t) interp1(t_data, u, t, 'linear', 'extrap');

[t,x]=ode45(@(t,x)odefun(t,x,u_input),[0 :0.0001 :20],[0;0]);

[phi11,~]=lsim(sys2,x(:,2),t);
[phi22,~]=lsim(sys2,x(:,1),t);
[phi33,~]=lsim(sys3,u,t);

y_dot= @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
phi1 = @(t) interp1(t_data, phi11, t, 'linear', 'extrap');
phi2 = @(t) interp1(t_data, phi22, t, 'linear', 'extrap');
phi3 = @(t) interp1(t_data, phi33, t, 'linear', 'extrap');

gamma=[0.5,0,0;0,0.5,0;0,0,0.5];

[t,theta]=ode45(@(t,theta)odefun2(t,theta,phi1,phi2,phi3,gamma,y_dot),[0:0.0001:20],[-0.5;0.3;0.2]);

m_hat=1./theta(:,3);
k_hat=m_hat.*theta(:,2);
b_hat=m_hat.*(theta(:,1)+am);

m=1.315*ones(length(t_data),1);
k=0.725*ones(length(t_data),1);
b=0.225*ones(length(t_data),1);

figure();
plot(t,m_hat,t,k_hat,t,b_hat);
hold on;
plot(t,m,'Color','b','LineStyle','--');
hold on;
plot(t,k,'Color','r','LineStyle','--');
hold on;
plot(t,b,'Color','y','LineStyle','--');
title('Εκτιμήσεις παραμέτρων');
legend({'$\hat{m}$','$\hat{k}$','$\hat{b}$','m','k','b'},'Interpreter','latex');
grid on;

m_hat = @(t) interp1(t_data, m_hat, t, 'linear', 'extrap');
b_hat = @(t) interp1(t_data, b_hat, t, 'linear', 'extrap');
k_hat = @(t) interp1(t_data, k_hat, t, 'linear', 'extrap');

x_hdot=theta(:,1).*phi11+theta(:,2).*phi22+theta(:,3).*phi33;
x_hat=lsim(sys4,x_hdot,t_data);

error=x(:,1)-x_hat;

figure();
plot(t,x(:,1),t,x_hat);
legend({'x(t)','$\hat{x}(t)$'},'Interpreter','latex');
grid on;

figure();
plot(t,error);
title('e=x(t)-$\hat{x}$(t)','Interpreter','latex');
grid on;

function dx=odefun(t,x,u)
    b=0.225;
    k=0.725;
    m=1.315;
    dx=[x(2);-(b/m)*x(2)-(k/m)*x(1)+(1/m)*u(t)];
end

function dtheta=odefun2(t,theta,phi1,phi2,phi3,gamma,y_dot)
    dtheta=gamma*(y_dot(t)-theta'*[phi1(t),phi2(t),phi3(t)]')*[phi1(t),phi2(t),phi3(t)]';
    disp(t);
end