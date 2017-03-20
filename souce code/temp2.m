[a, b] = textread('C:\Users\gzc\visual studio 2012\Projects\Project2\Project2\solution.txt');
plot(a,b,'o');
hold on
x=0:0.01:1;
plot(x, (exp(-5*(x-0.5).^2) - exp(-5/4))*(1 - exp(-5/4)))
legend('numerical', 'exact');
xlabel('x')
ylabel('Solution at y=0.5')
grid on