%% Satellite Orbit and Ground Track Simulation (LEO, MEO, GEO) - 3D animation
clc; clear; close all;
Re=6378.137; mu=398600.4418; we=7.2921159e-5;      % constants (km, s)
nm={'LEO','MEO','GEO'}; h=[700 20200 35786]; inc=[98 55 0]; c='rym';
a=Re+h; t=0:60:86164; th=we*t; nu=linspace(0,2*pi,361);
[x,y,z]=sphere(40);
g=@(la,lo,r) r*[cosd(la).*cosd(lo); cosd(la).*sind(lo); sind(la)];   % lat/lon -> xyz

%% Orbit coordinates, lat/lon, coverage
for k=1:3
    w=sqrt(mu/a(k)^3)*t; i=deg2rad(inc(k));
    P{k}=a(k)*[cos(w); sin(w)*cos(i); sin(w)*sin(i)];          % ECI
    ex= P{k}(1,:).*cos(th)+P{k}(2,:).*sin(th);                 % ECEF
    ey=-P{k}(1,:).*sin(th)+P{k}(2,:).*cos(th);
    la0{k}=asind(P{k}(3,:)/a(k)); lo0{k}=atan2d(ey,ex);        % sub-satellite point
    cov(k)=acosd(Re/a(k)*cosd(10))-10;                         % footprint angle (el = 10 deg)
    fprintf('%s: T = %.1f min, coverage radius = %.0f km\n',nm{k},2*pi*sqrt(a(k)^3/mu)/60,Re*deg2rad(cov(k)));
end

%% Figure 1: Earth + orbits + 3D orbit animation (inertial frame)
f1=figure; surf(Re*x,Re*y,Re*z,'FaceColor',[.2 .5 .9],'EdgeColor','none'); hold on;
axis equal; grid on; view(3); xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('3D Orbit Animation (Earth Model and Satellite Orbits)');
for k=1:3
    i=deg2rad(inc(k));
    hl(k)=plot3(a(k)*cos(nu),a(k)*sin(nu)*cos(i),a(k)*sin(nu)*sin(i),c(k),'LineWidth',1.5);
    hs(k)=plot3(NaN,NaN,NaN,'o','MarkerFaceColor',c(k),'MarkerEdgeColor','k','MarkerSize',8);
end
legend(hl,nm);

%% Figure 2: 3D ground track + coverage circle animation (on rotating-Earth frame)
f2=figure; az=linspace(0,360,181);
for k=1:3
    subplot(1,3,k); surf(Re*x,Re*y,Re*z,'FaceColor',[.2 .5 .9],'EdgeColor','none'); hold on;
    axis equal; grid on; view(3); axis(1.1*a(k)*[-1 1 -1 1 -1 1]);
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    title(sprintf('%s Ground Track + Coverage',nm{k}));
    hT(k)=plot3(NaN,NaN,NaN,c(k),'LineWidth',2); hC(k)=plot3(NaN,NaN,NaN,'k','LineWidth',1.5);
    hS(k)=plot3(NaN,NaN,NaN,'o','MarkerFaceColor',c(k),'MarkerEdgeColor','k','MarkerSize',8);
end
%% Endless animation of both figures (close either figure window to stop)
step=2; N=numel(t); j=1;            % smaller step / bigger pause = slower
while ishandle(f1) && ishandle(f2)
    for k=1:3
        set(hs(k),'XData',P{k}(1,j),'YData',P{k}(2,j),'ZData',P{k}(3,j));      % Figure 1 satellite
        p=la0{k}(j); q=lo0{k}(j); l=cov(k);
        ca=asind(sind(p)*cosd(l)+cosd(p)*sind(l)*cosd(az));                    % footprint circle
        co=q+atan2d(sind(az)*sind(l)*cosd(p),cosd(l)-sind(p)*sind(ca));
        T3=g(la0{k}(1:j),lo0{k}(1:j),Re*1.003);                                % track on globe
        C3=g(ca,co,Re*1.003); S3=g(p,q,a(k));                                  % circle, satellite
        set(hT(k),'XData',T3(1,:),'YData',T3(2,:),'ZData',T3(3,:));
        set(hC(k),'XData',C3(1,:),'YData',C3(2,:),'ZData',C3(3,:));
        set(hS(k),'XData',S3(1),'YData',S3(2),'ZData',S3(3));
    end
    drawnow limitrate; pause(0.03);
    j=mod(j-1+step,N)+1;            % wrap back to start -> endless loop
end