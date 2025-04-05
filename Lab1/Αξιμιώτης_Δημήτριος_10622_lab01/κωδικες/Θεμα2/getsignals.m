function [t,x,u]=getsignals()
% χρονικές εξισώσεις
[t,x]=ode45(@odefun,[0 :0.0001 :20],[0;0]);
u=4*cos(2*t);

function dx=odefun(t,x)

    % ορισμός παραμέτρων 
    m=0.75;
    L=1.25;
    c=0.15;
    g=9.81;
    A0=4;
    omega=2;

    % υλοποίηση διαφορικής εξίσωσης
    dx=[x(2);(-c/(m*L^2))*x(2)-(g/L)*x(1)+(1/(m*L^2))*A0*cos(omega*t)];

end

end