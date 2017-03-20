[x, p, u, rho] = textread('C:\Users\gzc\visual studio 2012\Projects\Project3\Project3\project3.txt');
figure(1)
plot(x,p,'o');
xlabel('position')
ylabel('Pressure')
grid on

figure(2)
plot(x,u,'o');
xlabel('position')
ylabel('Velocity')
grid on

figure(3)
plot(x,rho,'o');
xlabel('position')
ylabel('Density')
grid on
