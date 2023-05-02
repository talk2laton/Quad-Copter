function omegas2 = PDController(QuadCopter, Target)
    Ix = QuadCopter.MomIn(1); Iy = QuadCopter.MomIn(2); 
    Iz = QuadCopter.MomIn(3); l = QuadCopter.ArmLength; 
    k = QuadCopter.LiftCoefficient; b = QuadCopter.RotorDrag; 
    m = QuadCopter.Mass; s = QuadCopter.State; 
    kd = QuadCopter.AirResistance;
    K = [0.2     1.5
         0.2     1.5
         0.5     1.1
         0.2     2.6
         0.2     2.6
         0.5     2.6]; 
    P = K(:,1);  D = K(:,2); g = 9.81; [R, W] = Body_Intertia_Trans(s);
    sB = [zeros(3,1); zeros(3,1); R'*s(7:9); W*s(10:12)];
    
    if(numel(Target) == 3)
        T     = (m*g + P(3)*(Target(3) - s(3)) + (kd(3,3) - D(3))*s(9))/(cos(s(4))*cos(s(5))); %thrust
        theta = max(-0.5, min(0.5, (P(1)*(Target(1) - s(1)) + (kd(1,1)- D(1))*s(7))/T));
        phi   = max(-0.5, min(0.5, -(P(2)*(Target(2) - s(2)) + (kd(2,2) - D(2))*s(8))/T));
        xi    = 0;
    elseif(numel(Target) == 4)
        FI = [P(1)*(Target(1) - s(1)) + (kd(1,1) - D(1))*s(7);
              P(2)*(Target(2) - s(2)) + (kd(2,2) - D(2))*s(8);
              P(3)*(Target(3) - s(3)) + (kd(3,3) - D(3))*s(9) + m*g];
        FB = R'*FI; T = FI(3)/(cos(s(4))*cos(s(5)));
        Fx = FB(1)/T; Fy = FB(2)/T; Fz = FB(3)/T;
        xi = Target(4)-s(6); cxi = cos(xi); sxi = sin(xi); 
        phi   = max(-0.5, min(0.5, asin(Fx*sxi - Fy*cxi)));
        theta = max(-0.5, min(0.5, atan((Fx*cxi + Fy*sxi)/Fz)));
    end

    taux  = (P(4)*(phi - sB(4)) - D(4)*sB(10))*Ix; %roll
    tauy  = (P(5)*(theta - sB(5)) - D(5)*sB(11))*Iy; %pitch
    tauz  = (P(6)*(xi - sB(6)) - D(6)*sB(12))*Iz; % yaw
    omegas2 = [T/(4*k) + tauz/(4*b) + tauy/(4*k*l) + taux/(4*k*l);
               T/(4*k) - tauz/(4*b) + tauy/(4*k*l) - taux/(4*k*l);
               T/(4*k) + tauz/(4*b) - tauy/(4*k*l) - taux/(4*k*l);
               T/(4*k) - tauz/(4*b) - tauy/(4*k*l) + taux/(4*k*l)]; 
    omegas2 = max(0, omegas2)
end

function [R, W] = Body_Intertia_Trans(state)
    phi = state(4); theta = state(5); xi = state(6); 
    s1  = sin(xi); s2 = sin(theta); s3 = sin(phi); 
    c1  = cos(xi); c2 = cos(theta); c3 = cos(phi); 
    R   = [ c1*c2   c1*s2*s3-c3*s1   s1*s3+c1*c3*s2
            c2*s1   c1*c3+s1*s2*s3   c3*s1*s2-c1*s3
             -s2         c2*s3               c2*c3      ];

    W = [1      0     -s2
         0     c3    c2*s3
         0    -s3    c2*c3];

end


