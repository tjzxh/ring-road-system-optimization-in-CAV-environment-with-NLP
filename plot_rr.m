% ts = trackingScenario;
N = 22;T=100;Vlen=5;step=1;
height = 0;Rlen=230;r=Rlen/pi/2;
wayPoints = zeros(N,T,2,1);theta = zeros(N,T);
for t=2:T
    theta(1,t)=theta(1,t-1)+(sol.v(1,t-1)*step+sol.a(1,t-1)*step^2/2+Vlen)/r;
end
for n=2:N
    for tt=1:T
        theta(n,tt)=theta(n-1,tt)-(sol.s(n-1,tt)+Vlen)/r;
    end
end
wayPoints(1:N,1:T,1,:) = r*cos(theta(1:N,1:T));
wayPoints(1:N,1:T,2,:) = r*sin(theta(1:N,1:T));
sz=100;
% scatter(wayPoints(:,:,1,1),wayPoints(:,:,2,1))
plot((r-5)*cos(0:pi/50:2*pi),(r-5)*sin(0:pi/50:2*pi),(r+5)*cos(0:pi/50:2*pi),(r+5)*sin(0:pi/50:2*pi),Color='k',LineWidth=2);
pbaspect([1 1 1])
daspect([1 1 1])
% n1=1:7;n2=9:N-1;HV=[n1,n2];AV=[8,N];
% HV=[1,2,4:7,9:N-1];AV=[3,8,N];
% HV=[1,2:3,5:8,11:14,16:N-1];AV=[4,9,10,15,N];
HV=[1,4,6,7,11:12,14:16,19:N-1];AV=[2,3,5,8,9,10,13,17,18,N];
hold on
for t = 1:T
    h = scatter(wayPoints(HV,t,1,1),wayPoints(HV,t,2,1),sz,'o','filled','black');            % draw something on the trajectory
    targ = scatter(wayPoints(AV,t,1,1),wayPoints(AV,t,2,1),sz,'o','filled','red');
    exportgraphics(gcf,'testAnimated.gif','Append',true);
%     pause(0.05)                                % wait a minute
    delete(h)                                 % delete it
    delete(targ)
end
hold off

% x = 0:0.01:1;
% p = plot(nan,nan);
% p.XData = x;
% for n = 1:0.5:5
%       p.YData = x.^n;
%       exportgraphics(gcf,'testAnimated.gif','Append',true);
% end