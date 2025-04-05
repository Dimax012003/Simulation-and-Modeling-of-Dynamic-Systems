A=[0.01,0.1,0.2,0.5,2.5,5,10,50,100];

L_bar=[];
m_bar=[];
c_bar=[];
for i=1:length(A)
    theta1=create_error_graph_7(A(i));
    L_bar=[L_bar theta1(1)];
    m_bar=[m_bar theta1(2)];
    c_bar=[c_bar theta1(3)];
    disp(i);
end

figure();
plot(A,L_bar,LineWidth=1);
hold on;
plot(A,m_bar,LineWidth=1);
hold on;
plot(A,c_bar,LineWidth=1);
hold off;
xlabel('A');
ylabel('Estimation error');
legend({'$\tilde{L}$','$\tilde{m}$','$\tilde{c}$'},'Interpreter','latex');
grid on;

function [theta]=create_error_graph_7(A)

[t,x,u]=getsignals_A(A);

sys1=tf([-1 0],[1 3 2]);
sys2=tf([0 -1],[1 3 2]);
sys3=tf([0 1],[1 3 2]);

Ts=0.1;
Dt=0.0001;
step=Ts/Dt;

k=[1:step:length(t)]';
t_new=[0:Ts:20]';

[x1,~]=lsim(sys1,x(k,1),t_new);
[x2,~]=lsim(sys2,x(k,1),t_new);
[x3,~]=lsim(sys3,u(k),t_new);

X=[x1,x2,x3];
y=x(k,1);

[theta,~]=least_squares(X,y);

L_hat = 9.81 / (theta(2)+2);
m_hat = 1 / (theta(3)*L_hat^2);
c_hat = (3+theta(1))*(m_hat*L_hat^2);



L_bar=1.25-L_hat;
m_bar=0.75-m_hat;
c_bar=0.15-c_hat;


theta=[L_bar,m_bar,c_bar];
disp(theta);

function dq=odefun(t,q,m_hat,L_hat,c_hat,A0)

    % ορισμός παραμέτρων 
    g=9.81;
    omega=2;
    % υλοποίηση διαφορικής εξίσωσης
    dq=[q(2);(-c_hat/(m_hat*L_hat^2))*q(2)-(g/L_hat)*q(1)+(1/(m_hat*L_hat^2))*A0*cos(omega*t)];

end
end