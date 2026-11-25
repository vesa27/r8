clc; clear; close all;
Re=6378.137; mu=398600.4418; we=7.2921159e-5;
nm={'LEO','MEO','GEO'}; h=[700 20200 35786]; inc=[98 55 0]; c='rym';
a=Re+h; t=0:60:86164; th=we*t; nu=linspace(0,2*pi,361);
[x,y,z]=sphere(40);

for k=1:3
    w=sqrt(mu/a(k)^3)*t; i=deg2rad(inc(k));
    P{k}=a(k)*[cos(w); sin(w)*cos(i); sin(w)*sin(i)];
    ex= P{k}(1,:).*cos(th)+P{k}(2,:).*sin(th);
    ey=-P{k}(1,:).*sin(th)+P{k}(2,:).*cos(th);
    la0{k}=asind(P{k}(3,:)/a(k)); lo0{k}=atan2d(ey,ex);
    cov(k)=acosd(Re/a(k)*cosd(10))-10;
    fprintf('%s: T = %.1f min, coverage radius = %.0f km\n',nm{k},2*pi*sqrt(a(k)^3/mu)/60,Re*deg2rad(cov(k)));
end

figure('Position',[100 50 900 850]); az=linspace(0,360,181);
for k=1:3
    lo=lo0{k}; lo([false abs(diff(lo))>180])=NaN;
    p=la0{k}(1); q=lo0{k}(1); l=cov(k);
    ca=asind(sind(p)*cosd(l)+cosd(p)*sind(l)*cosd(az));
    co=mod(q+atan2d(sind(az)*sind(l)*cosd(p),cosd(l)-sind(p)*sind(ca))+180,360)-180;
    co([false abs(diff(co))>180])=NaN;
    subplot(3,1,k);
    plot(lo,la0{k},c(k),co,ca,'k',q,p,'ko','MarkerFaceColor',c(k)); grid on;
    axis([-180 180 -90 90]); xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
    title(sprintf('%s Ground Track + Coverage',nm{k}));
end

f1=figure; surf(Re*x,Re*y,Re*z,'FaceColor',[.2 .5 .9],'EdgeColor','none'); hold on;
axis equal; grid on; view(3); xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('3D Orbit Animation (Earth Model and Satellite Orbits)');
for k=1:3
    i=deg2rad(inc(k));
    hl(k)=plot3(a(k)*cos(nu),a(k)*sin(nu)*cos(i),a(k)*sin(nu)*sin(i),c(k),'LineWidth',1.5);
    hs(k)=plot3(NaN,NaN,NaN,'o','MarkerFaceColor',c(k),'MarkerEdgeColor','k','MarkerSize',8);
end
legend(hl,nm);

step=2; N=numel(t); j=1;
while ishandle(f1)
    for k=1:3
        set(hs(k),'XData',P{k}(1,j),'YData',P{k}(2,j),'ZData',P{k}(3,j));
    end
    drawnow limitrate; pause(0.03);
    j=mod(j-1+step,N)+1;
end
