close all
camva(5); camlight(30,30); 
fig = gcf; fig.Color = 'w'; 
ax = gca;
Q = QuadCopter(ax, [-15,  15,  0, 0], 0.2);
TargetList = [-15,   15,   10,  pi/2
               15,   15,   10,  pi/2
               15,  -15,   10,  0
              -15,  -15,   10,  -pi/2
              -15,  -15,    0,  -pi/2];
axis([-20,20,-20,20,-2,14]); grid on; hold on
FlightSimulator(Q, TargetList, 30, 0.2, 'flightcontrol11.avi');
