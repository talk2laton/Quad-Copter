close all
camva(5); camlight(30,30); 
fig = gcf; fig.Color = 'w'; 
fig.Position = [962 42 958 954];
ax = gca;
Q = QuadCopter(ax = ax, pos = [-15,  15,  0, pi/2], scalefactor= 0.2);
TargetList = [-15,   15,   10, 0
               15,   15,   10, -pi/2
               15,  -15,   10, -pi/2
              -15,  -15,   10, 0
              -15,  -15,    0, pi/2];
axis([-25,25,-25,25,-2,14]); grid on; hold on
FlightSimulator(Q, TargetList, 30, 0.01, 'flightsimulation2.avi') ; 
