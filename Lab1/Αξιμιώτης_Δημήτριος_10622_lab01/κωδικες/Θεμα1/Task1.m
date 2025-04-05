[t,x,u]=getsignals();
% plots
figure();
plot(t,x(:,1),t,x(:,2),t,u);
xlabel('t(s)');
ylabel('State Variables and input signal');
legend('q(t)','dq(t)/dt','u');
grid on;