%class definition
classdef QuadCopter < handle
    properties (SetAccess = immutable)
        ArmLength       = 0.225;
        LiftCoefficient = 2.980e-5;
        RotorDrag       = 1.140e-7;
        Mass            = 0.468;
        MomIn           = [4.856e-3;4.856e-3;8.801e-3];
        GyroConstant    = 3.357e-5;
        AirResistance   = 0.25*eye(3);
    end

    properties(SetAccess = private)
        FlightTime  (1,1)  double
        State       (12,1) double
        GPS         (3,1)  double
        GYRO        (3,1)  double
        Update      
        Dynamics    
    end

    properties (Access = private)
        TotalRotorAngle (4,1)  double = zeros(4,1);
        Axis            (4,4)  double = eye(4); 
        PosHg
        RotHg
        Spin1
        Spin2
        Spin3
        Spin4
    end

    methods 
        function Q = QuadCopter(options)
            arguments
                options.ax  (1,1) matlab.graphics.axis.Axes
                options.pos (1,4) double
                options.scalefactor (1,1) double
                options.parametersfile (1,1) string
                options.userdynamics (1,1) function_handle
            end
            [ax, pos, scalefactor] = setArguments(options);
            axis equal;
            C = colororder; 
            camlight(ax, 'left', 'infinite ');
            camlight(ax, 'right', 'infinite ');
            camlight(ax, 'headlight', 'infinite ');
            Q.PosHg = hgtransform('Parent',ax); 
            set(Q.PosHg, 'Matrix', makehgtform('translate', pos(1:3)));
            Size = hgtransform('Parent', Q.PosHg); 
            set(Size, 'Matrix', makehgtform('scale', scalefactor));
            Q.RotHg = hgtransform('Parent', Size); 
            Q.State = [pos(1:3)';zeros(9,1)]; Q.State(6) = pos(4);
            Q.GPS   = pos(1:3);  % center of gravity of the helicopter
            Q.GYRO  = [0,0,pos(4)]; % be replaced with get GYRO orientation
            set(Q.RotHg, 'Matrix', makehgtform('zrotate', pos(4)));
            Q.FlightTime = 0; 
            c12 = C(1:2, :); c34 = C(3:4, :); c56 = C(5:6, :); c7 = C(7, :);
            [hing1, hing2, hing3, hing4] = MakeBody(Q.RotHg, c12); 
            Q.Spin1 = MakeElectricMotor(MakeArm(hing1, c34), c56); 
            Q.Spin2 = MakeElectricMotor(MakeArm(hing2, c34), c56);
            Q.Spin3 = MakeElectricMotor(MakeArm(hing3, c34), c56);
            Q.Spin4 = MakeElectricMotor(MakeArm(hing4, c34), c56);
            MakeRotors(Q.Spin1, -1, c7); MakeRotors(Q.Spin2, 1, c7);
            MakeRotors(Q.Spin3, -1, c7); MakeRotors(Q.Spin4, 1, c7);
            Q.Update   = @(wsqr, ds, dt) update(Q, wsqr, ds, dt);
            if(isfield(options, 'userdynamics') && ...
                    isa(options.userdynamics,'function_handle'))
                Q.Dynamics = @(wsqr, dt) options.userdynamics(Q, wsqr, dt);
            else
                Q.Dynamics = @(wsqr, dt) dynamics(Q, wsqr, dt);
            end
            if(isfield(options, 'parametersfile'))
                param = readtable(options.parametersfile);
                param.Properties.RowNames = param{:, "Var1"}; 
                param(:, "Var1") = [];
                Q.ArmLength       = param{"ArmLength", 1};
                Q.LiftCoefficient = param{"LiftCoefficient", 1};
                Q.RotorDrag       = param{"RotorDrag", 1};
                Q.Mass            = param{"Mass", 1};
                Q.MomIn           = param{"MomIn", 1:3};
                Q.GyroConstant    = param{"GyroConstant", 1};
                Q.AirResistance   = diag(param{"AirResistance", 1:3});
            end
        end
        
        function ds = dynamics(Q, wsqr, dt)
            k = Q.LiftCoefficient; w = sqrt(wsqr);
            l = Q.ArmLength; b = Q.RotorDrag; J = Q.GyroConstant;
            kd = Q.AirResistance; m = Q.Mass; Ixx   = Q.MomIn(1) ;  
            Iyy = Q.MomIn(2);  Izz = Q.MomIn(3);  g = 9.81;
            T    = k * (wsqr(1) + wsqr(2) + wsqr(3) + wsqr(4)); 
            taux = k*l*(wsqr(1) - wsqr(2) - wsqr(3) + wsqr(4));
            tauy = k*l*(wsqr(1) + wsqr(2) - wsqr(3) - wsqr(4));
            tauz = b * (wsqr(1) - wsqr(2) + wsqr(3) - wsqr(4));
            wr = -w(1) + w(2) - w(3) + w(4);
            function ds = dsdt(s)
                vI  = s(7:9); wI  = s(10:12); 
                [R, W, invW, dinvW] = Body_Intertia_Trans(s);
                wB = W*wI;
                p    = wB(1); q = wB(2);   r  = wB(3); 
                vI_p = [0;0;-g] + (R*[0;0;T] - kd*vI)/m;
                wB_p = [((Iyy - Izz)*(q*r) - J*q*wr + taux)/Ixx
                        ((Izz - Ixx)*(p*r) + J*p*wr + tauy)/Iyy
                        ((Ixx - Iyy)*(p*q) + tauz)/Izz];
                wI_p = dinvW*wB + invW*wB_p;
                ds   = [vI; wI; vI_p; wI_p];
            end
            s = Q.State; 
            K1 = dt*dsdt(s); K2 = dt*dsdt(s + 0.5*K1); 
            K3 = dt*dsdt(s + 0.5*K2); K4 = dt*dsdt(s + K3); 
            ds = (K1 + 2*K2 + 2*K3 + K4)/6;
        end
        
        function update(Q, wsqr, ds, dt) 
            Q.State = Q.State + ds;
            Q.GPS = Q.State(1:3); 
            Q.GYRO = Q.State(4:6); 
            Q.FlightTime = Q.FlightTime + dt;
            set(Q.PosHg,'Matrix', makehgtform('translate',Q.GPS));
            Q.Axis = makehgtform('zrotate', Q.GYRO(3), 'yrotate', Q.GYRO(2), 'xrotate', Q.GYRO(1));
            set(Q.RotHg,'Matrix', Q.Axis);
            RotorsRotate(Q, sqrt(wsqr), dt); 
            drawnow
        end
    end
end

function [ax, pos, scalefactor] = setArguments(options)
    if isfield(options, 'scalefactor')
        scalefactor = options.scalefactor;
    else
        scalefactor = 1;
    end
    if isfield(options, 'pos')
        pos = options.pos;
    else
        pos = [0,0,4,0];
    end
    if isfield(options, 'ax')
        ax = options.ax;
    else
        figure(Color = 'w');
        ax = gca; daspect([1,1,1]);
    end
end

function [hing1, hing2, hing3, hing4] = MakeBody(hg, clr)
    t = linspace(0,2*pi,51);
    x = sign(cos(t)).*abs(cos(t)).^0.2;
    y = sign(sin(t)).*abs(sin(t)).^0.2;
    X = [4;4;2;2;4;4]*x; Y = [6;6;4;4;6;6]*y; 
    Z = [3;2;2;-2;-2;-3]*ones(size(t));
    surf(X, Y, Z, 'FaceColor', clr(1,:), 'EdgeAlpha', 0.1, Parent = hg); 
    [X,Y] = meshgrid([-5,-4,4,5],[-7,7]);
    Z = 0*X; Z(X==-5) = 2.5; Z(X==-4) = 3; Z(X==5) = 2.5; Z(X==4) = 3;
    surf(X, Y, Z, 'FaceColor', clr(2,:), 'EdgeAlpha', 0.1, Parent = hg); 
    surf(X, Y, Z+0.2, 'FaceColor', clr(2,:), 'EdgeAlpha', 0.1, Parent = hg); 
    surf(X, Y, -Z, 'FaceColor', clr(2,:), 'EdgeAlpha', 0.1, Parent = hg); 
    surf(X, Y, -Z-0.2, 'FaceColor', clr(2,:), 'EdgeAlpha', 0.1, Parent = hg); 
    x = [X(1,:), flip(X(1,:))]; y = [Y(1,:), flip(Y(1,:))]; 
    z = [Z(1,:), flip(Z(1,:)+0.2)]; 
    fill3(x, y, z, clr(2,:), Parent = hg); fill3(x, y+14, z, clr(2,:), Parent = hg);
    fill3(x, y, -z, clr(2,:), Parent = hg); fill3(x, y+14, -z, clr(2,:), Parent = hg);
    x = [X(:,end); flip(X(:,end))]; y = [Y(:,end); flip(Y(:,end))]; 
    z = [Z(:,end); flip(Z(:,end)+0.2)]; 
    fill3(x, y, z, clr(2,:), Parent = hg); fill3(x-10, y, z, clr(2,:), Parent = hg);
    fill3(x, y, -z, clr(2,:), Parent = hg); fill3(x-10, y, -z, clr(2,:), Parent = hg);
    FixPos1 = hgtransform('Parent',hg); set(FixPos1,'Matrix', makehgtform('translate',[3,5,0]));
    FixPos2 = hgtransform('Parent',hg); set(FixPos2,'Matrix', makehgtform('translate',[3,-5,0]));
    FixPos3 = hgtransform('Parent',hg); set(FixPos3,'Matrix', makehgtform('translate',[-3,5,0]));
    FixPos4 = hgtransform('Parent',hg); set(FixPos4,'Matrix', makehgtform('translate',[-3,-5,0]));
    FixRot1 = hgtransform('Parent',FixPos1); set(FixRot1,'Matrix', makehgtform('zrotate',atan2(5,3)));
    FixRot2 = hgtransform('Parent',FixPos2); set(FixRot2,'Matrix', makehgtform('zrotate',atan2(-5,3))); 
    FixRot3 = hgtransform('Parent',FixPos4); set(FixRot3,'Matrix', makehgtform('zrotate',atan2(-5,-3)));
    FixRot4 = hgtransform('Parent',FixPos3); set(FixRot4,'Matrix', makehgtform('zrotate',atan2(5,-3)));  
    hing1 = hgtransform('Parent',FixRot1); hing2 = hgtransform('Parent',FixRot2); 
    hing3 = hgtransform('Parent',FixRot3); hing4 = hgtransform('Parent',FixRot4); 
end
   
function Fix = MakeArm(hg, clr)
    [Xc, Yc, Zc] = cylinder([0.5, 1, 1, 0.5]); Zc(1:2,:) = 0; Zc(3:4,:) = 4;
    N = 11; t = linspace(0,pi,N);
    y = sign(cos(t)).*abs(cos(t)).^0.2;
    z = 4*sign(sin(t)).*abs(sin(t)).^0.2;
    t = flip(t);
    y = [y, 0.8*sign(cos(t)).*abs(cos(t)).^0.2, y(1)];  
    z = [z, 3.5*sign(sin(t)).*abs(sin(t)).^0.2, z(1)]; 
    x = 0.1*(0:100)'; h = (2.7*x).^(1/3);
    X = x*ones(1,2*N+1); Y = repmat(y, size(x)); 
    Z = max(repmat(z, size(x)), h*ones(1,2*N+1));
    surf(Xc, Yc, Zc-2, 'FaceColor', clr(1,:), 'EdgeAlpha', 0.1, Parent = hg); 
    Zc(1:2,:) = 2; 
    surf(10+Xc, Yc, Zc-2, 'FaceColor', clr(1,:), 'EdgeAlpha', 0.1, Parent = hg);
    surf(X, Y, Z-2, 'FaceColor', clr(1,:), 'EdgeAlpha', 0.1, Parent = hg); 

    xyz = load('xyz.txt');
    L = cumsum([0; vecnorm(diff(xyz),2,2)]); L = L/L(end)*1000;
    xyz(diff(L) == 0,:)=[]; L(diff(L) == 0,:) = [];
    [X,Y,Z] = tube(@(t)0.3+0*t,20,@(t)interp1(L, xyz, t),[0,1000],201);
    [X,Y,Z] = rotate(X, Y, 0.5+Z, [0,0,1],-pi/4);
    surf(10+X, Y, Z, 'FaceColor', 'w', 'EdgeAlpha', 0.1, Parent = hg); 
    [X,Y,Z] = tube(@(t)0.3+0*t,20,@(t)[8*cos(t), 8*sin(t), 0*t],[-pi/4,pi/4],31);
    surf(10+X, Y, Z, 'FaceColor', 'w', 'EdgeAlpha', 0.1, Parent = hg); 
    surf(10+X, Y, Z-2, 'FaceColor', 'w', 'EdgeAlpha', 0.1, Parent = hg); 
    surf(10+X, Y, Z+2, 'FaceColor', 'w', 'EdgeAlpha', 0.1, Parent = hg); 

    Fix = hgtransform('Parent',hg); 
    set(Fix,'Matrix', makehgtform('translate',[10,0,2.4]));
end

function Spin = MakeElectricMotor(hg, clr)
    [Xc, Yc, Zc] = cylinder([0, 0.7, 0.7]); Zc(1:2,:) = -0.2; Zc(3,:) = -2;
    surf(Xc, Yc, Zc, 'FaceColor', clr(1,:), 'EdgeAlpha', 0.1, Parent = hg); 
    surf(0.2*Xc, 0.2*Yc, Zc+0.3, 'FaceColor', clr(2,:), 'EdgeAlpha', 0.1, Parent = hg); 
    Spin = hgtransform('Parent',hg); 
end

function MakeRotors(hg, sgn, clr)
    t = linspace(0,2*pi,31); c = cos(t); ac = abs(c); s = sin(t); as = abs(s);
    B = 2; T = 0.2; C = 0.05; P = 1; E = 1; y2 = 0.5+0.5*ac.^B./c-0.08;
    z2 = T/2*as.^B./s.*(1-y2.^P) + C.*sin(y2.^E*pi); z2(1) = z2(end);
    t = atan2(z2, y2); c = cos(t); s = sin(t);
    y1 = 0.05*c; z1 = 0.05*s; x1 = ones(size(t));
    X = [0;0.7;1;6]*x1; Y = sgn*[y1;y1;y2;y2]; Z = [z1;z1;z2;z2];
    for n = 1:3
        [X, Y, Z] = rotate(X, Y, Z, [0,0,1], 2*pi/3);
        surf(X, Y, Z, 'FaceColor', clr, 'EdgeAlpha', 0.1, Parent = hg); 
        fill3(X(end,:), Y(end,:), Z(end,:), clr, Parent = hg); 
    end
end

function RotorsRotate(Q, Omegas, dt)
    Q.TotalRotorAngle = Q.TotalRotorAngle + Omegas*dt;
    set(Q.Spin1, 'Matrix', makehgtform('zrotate', Q.TotalRotorAngle(1)));
    set(Q.Spin2, 'Matrix', makehgtform('zrotate', -Q.TotalRotorAngle(2)));
    set(Q.Spin3, 'Matrix', makehgtform('zrotate', Q.TotalRotorAngle(3)));
    set(Q.Spin4, 'Matrix', makehgtform('zrotate', -Q.TotalRotorAngle(4)));
end 

function [X, Y, Z] = rotate(X, Y, Z, U, theta)
    M = Rxyz(U,theta);
    for i = 1:size(X,1)
        xyz= M*[X(i,:);Y(i,:);Z(i,:)];
        X(i,:) = xyz(1,:);
        Y(i,:) = xyz(2,:);
        Z(i,:) = xyz(3,:);
    end
end

function [X,Y,Z] = tube(Rfxn,Nr,Paramtricfxn,tbound,Nt)
    T = linspace(tbound(1),tbound(2),Nt + 1)'; Trans = Paramtricfxn(T);
    R = Rfxn(T); [X,Y,Z] = cylinder(R,Nr); Z = 0*Z;
    [axldir, angle] = PathDir(Paramtricfxn,T); 
    for n = 1:Nt + 1
        R = Rxyz(axldir(n,:),angle(n));
        mat = R*[X(n,:);Y(n,:);Z(n,:)];
        X(n,:) = Trans(n,1) + mat(1,:);
        Y(n,:) = Trans(n,2) + mat(2,:);
        Z(n,:) = Trans(n,3) + mat(3,:);
    end
end

function [dir, ang] = PathDir(Paramtricfxn,T)
    T(end) = T(end) -1e-6; Tpdt = T + 1e-6;
    del  = diff(cat(3, Paramtricfxn(T), Paramtricfxn(Tpdt)),1,3);
    v    = vecnorm(del, 2, 2); 
    V    = del./(v*ones(1,3));
    dir  = [- V(:,2), V(:,1), zeros(size(V,1),1)];
    for j = find(vecnorm(dir, 2, 2) == 0)'
        dir(j,:) = dir(j-1,:);
    end
    v    = vecnorm(dir, 2, 2); 
    dir  = dir./(v*ones(1,3));
    ang  = acos(V(:,3));
end

function R = Rxyz(U,theta)
    c = cos(theta); s = sin(theta); ux = U(1); uy = U(2); uz = U(3);
    R = [ux*ux*(1 - c) +    c  uy*ux*(1 - c) - uz*s  uz*ux*(1 - c) + uy*s
         uy*ux*(1 - c) + uz*s  uy*uy*(1 - c) +    c  uz*uy*(1 - c) - ux*s
         uz*ux*(1 - c) - uy*s  uy*uz*(1 - c) + ux*s  uz*uz*(1 - c) +    c];  
end

function [R, W, invW, dinvW] = Body_Intertia_Trans(state)
    phi = state(4); theta = state(5); xi = state(6); 
    phip = state(10); thetap = state(11); 
    s1 = sin(xi); s2 = sin(theta); s3 = sin(phi); 
    c1 = cos(xi); c2 = cos(theta); c3 = cos(phi); 
    t2 = tan(theta);
    R = [c1*c2   c1*s2*s3-c3*s1   s1*s3+c1*c3*s2
         c2*s1   c1*c3+s1*s2*s3   c3*s1*s2-c1*s3
          -s2        c2*s3            c2*c3     ];

    W    = [1       0       -s2
            0      c3      c2*s3
            0     -s3      c2*c3];

    invW = [1     s3*t2    c3*t2
            0      c3       -s3
            0     s3/c2    c3/c2];

    dinvW = [0   phip*c3*t2 + thetap*s3/c2^2   -phip*s3*t2 + thetap*c3/c2^2
             0           -phip*s3                       -phip*c3
             0   phip*c3/c2 + thetap*s3*t2/c2  -phip*s3/c2 + thetap*c3*t2/c2];
end