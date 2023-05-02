function FlightSimulator(QuadCopter, TargetList, framerate, Tol, videoname)
    clc;
    global vidObj
    writevideo = exist('videoname','var');
    view(-42.8833, 7.8978); t = 0; dt = 1/framerate; 
    %% Start count down
    plot3(TargetList(:,1), TargetList(:,2), TargetList(:,3), 'or', 'MarkerFaceColor','r');
    for n = 1:size(TargetList,1)
        Target = TargetList(n,:);
        text(Target(1),Target(2),Target(3)+5,num2str(n), "FontSize",12);
    end
    CountDownFrom(5);
    Path = [];
    %% Capturing Video
    if(writevideo)
        vidObj = VideoWriter(videoname); set(vidObj, 'FrameRate',framerate); open(vidObj);
    end
    Totaltime = 0;
    %% Navigating
    errorGYRO = 0;
    for n = 1:size(TargetList,1)
        Target = TargetList(n,:)'; 
        errorGPS = norm(Target(1:3) - QuadCopter.State(1:3));
        if(numel(Target) == 4)
            errorGYRO = abs(Target(4) - QuadCopter.State(6));
        end
        iter = 0; 
        Path = [Path, QuadCopter.GPS];
        while ((errorGPS > Tol || errorGYRO > Tol/10) && iter < 5000)
            wsqr = PDController2(QuadCopter, Target);
            ds = QuadCopter.Dynamics(wsqr, dt); 
            QuadCopter.Update(wsqr, ds, dt);
            Path = [Path, QuadCopter.GPS];
            r = QuadCopter.State(1:3); v = QuadCopter.State(7:9);
            title({sprintf('$Going~to~point (n = %2d) = [%8.2f,~ %8.2f,~ %8.2f]$,',n, Target(1), Target(2), Target(3)),...
                   sprintf('$Current~Position: r (m) = [%8.2f,~ %8.2f,~ %8.2f]$,',r(1), r(2), r(3)),...
                   sprintf('$Current~Velocity: \\dot{r} (m/s) = [%8.2f,~ %8.2f,~ %8.2f]$, ',v(1), v(2), v(3))}, ...
                   'interpreter','latex');
            drawnow; iter = iter + 1; 
            errorGPS = norm(Target(1:3) - QuadCopter.State(1:3));
            if(numel(Target) == 4)
                errorGYRO = abs(Target(4) - QuadCopter.State(6));
            end
            if(writevideo)
                image = getframe(gcf);
                writeVideo(vidObj,image);
            end
            t = t + dt;
        end
        
        %% Hover mode
        disp('===============================================================');
        disp(QuadCopter.FlightTime);
        disp('===============================================================');
    end
    plot3(Path(1,:), Path(2,:), Path(3,:))
    if(writevideo)
        close(vidObj);
    end
end

function CountDownFrom(N)
    i = N;
    for n = 1:N
        pause(1);
        title(num2str(i));
        i = i - 1;
    end
    pause(1); title('Start');
end

