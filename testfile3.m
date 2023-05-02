close all
camva(10); camlight(30,30); 
fig = gcf; fig.Color = 'w'; 
ax = gca;
Q = QuadCopter(ax, [0,  0,  0, 0], 0.4);
TargetList = [-30+60*rand(10,2), 100*rand(10,1), pi*(rand(10,1)-0.5)];
axis([-40,40,-40,40,-10,110]); grid on; hold on
FlightSimulator(Q, TargetList, 30, 0.5, 'flightcontrol15.avi') ; 