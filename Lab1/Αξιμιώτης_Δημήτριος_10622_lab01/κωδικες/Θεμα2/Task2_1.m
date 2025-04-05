e=create_error_graph_1(0.1);
function e=create_error_graph_1(Ts)
[t,x,u]=getsignals();

% plots
figure();
plot(t,x(:,1),t,x(:,2));
xlabel('t(s)');
ylabel('State Variables');
legend('q(t)','dq(t)/dt');
grid on;

sys1=tf([-1 0],[1 1]);
sys2=tf([0 -1],[1 1]);
sys3=tf([0 1],[1 1]);

Dt=0.0001;
step=Ts/Dt;

k=[1:step:length(t)]';
t_new=[0:Ts:20]';

[x1,~]=lsim(sys2,x(k,2),t_new);
[x2,~]=lsim(sys2,x(k,1),t_new);
[x3,~]=lsim(sys3,u(k),t_new);

X=[x1,x2,x3];
y=x(k,2);

[theta,~]=least_squares(X,y);

L_hat = 9.81 / theta(2);
m_hat = 1 / (theta(3)*L_hat^2);
c_hat = (1+theta(1))*(m_hat*L_hat^2);

[t,q_hat]=ode45(@(t,q_hat)odefun(t,q_hat,m_hat,L_hat,c_hat),[0 :Dt :20],[0;0]);
disp("    L         m         c");
disp([L_hat,m_hat,c_hat]);

e=x(:,1)-q_hat(:,1);
figure();
plot(t,e);
title("Error plot");
xlabel('t');
ylabel('e');
grid on;

figure();
plot(t,x(:,1),t,q_hat(:,1),LineWidth=0.8,LineStyle="-");
xlabel('t');
ylabel('Angular values');
legend('Real angle','Estimated angle');
grid on;

figure();
plot(t,x(:,2),t,q_hat(:,2),LineWidth=0.8,LineStyle="-");
xlabel('t');
ylabel('Velocity values');
legend('Real velocity','Estimated velocity');
grid on;
disp(['mse=',num2str(mean(e.^(2)))]);
function dq=odefun(t,q,m_hat,L_hat,c_hat)

    % ορισμός παραμέτρων 
    g=9.81;
    A0=4;
    omega=2;
    % υλοποίηση διαφορικής εξίσωσης
    dq=[q(2);(-c_hat/(m_hat*L_hat^2))*q(2)-(g/L_hat)*q(1)+(1/(m_hat*L_hat^2))*A0*cos(omega*t)];

end
end