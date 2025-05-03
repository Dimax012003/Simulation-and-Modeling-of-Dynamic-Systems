[t,x]=ode45(@odefun,[0 :0.0001 :20],[0;0]);

t_data=t;
u=2.5*sin(t_data);

y_dot= @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
y = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');
u = @(t) interp1(t_data, u, t, 'linear', 'extrap');

gamma=[10,10,10];

[t,Z]=ode45(@(t,Z)odefun2(t,Z,u,y,y_dot,gamma,2),[0 :0.0001 :20],[0;0;0.2;0.2;0.2]);

theta1=Z(:,3);
theta2=Z(:,4);
theta3=Z(:,5);

m_hat=1./theta3;
k_hat=m_hat.*theta2;
b_hat=m_hat.*theta1;

x_hat1=Z(:,1);
x_hat2=Z(:,2);

figure();
plot(t,m_hat,t,k_hat,t,b_hat);
legend('m_hat','k_hat','b_hat');
grid on;

figure();
plot(t,x(:,1),t,x_hat1);
legend('x(t)','hat{x(t)}');
grid on;

error=x(:,1)-x_hat1;

figure();
plot(t,error);
legend('e');
grid on;

function dx=odefun(t,x)

    % ορισμός παραμέτρων 
    b=0.225;
    k=0.725;
    m=1.315;
    % υλοποίηση διαφορικής εξίσωσης
    dx=[x(2);-(b/m)*x(2)-(k/m)*x(1)+(1/m)*2.5*sin(t)];

end

function dZ=odefun2(t,Z,u,y,y_dot,gamma,n0)
    x = Z(1:2);
    theta = Z(3:5);
    n=n0*sin(2*pi*20);
    dx=[x(2);-theta(1)*x(2)-theta(2)*x(1)+theta(3)*2.5*sin(t)];
    dtheta(1)=-gamma(1)*(y_dot(t)-x(2))*x(2);
    dtheta(2)=-gamma(2)*(y_dot(t)-x(2))*x(1);
    dtheta(3)=gamma(3)*(y_dot(t)-x(2))*u(t);
    dZ=[dx;dtheta'];
    disp(t);
end