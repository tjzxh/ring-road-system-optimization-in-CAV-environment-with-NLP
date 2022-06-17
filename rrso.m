%% Define const
function ringroad=rrso(init_a,init_v,init_s,result)
N = 22;T=100;ts=1;AVN=10;
HV=[4,6,7,11:12,14:16,19:N-1];AV=[2,3,5,8,9,10,13,17,18,N];
Rlen=230;Vlen=5;
a = optimvar('a',N,T,'LowerBound',-5,'UpperBound',5);
v = optimvar('v',N,T,'LowerBound',0,'UpperBound',120/3.6);
s = optimvar('s',N,T,'LowerBound',1,'UpperBound',Rlen-N*Vlen);
e1 = optimvar('e1',N,T,'LowerBound',0,'UpperBound',2);
e2 = optimvar('e2',N,T,'LowerBound',0,'UpperBound',1);
ringroad = optimproblem('ObjectiveSense','maximize');
%% IDM
alpha=0.73;beta=1.67;v0=30;hd=1.5;s0=2;t=1:T;
%     aconsup = optimconstr(N-1,T-1);
%     aconsup(n,t)=a(n,t)<=alpha*(1-(v(n,t)./v0).^4-((s0+v(n,t).*hd+v(n,t).*(v(n,t)-v(n-1,t))/2/(alpha*beta).^0.5)./s(n,t)).^2)+e;
%     aconsup(1,t)=a(1,t)<=alpha*(1-(v(1,t)./v0).^4-((s0+v(1,t).*hd+v(1,t).*(v(1,t)-v(N,t))/2/(alpha*beta).^0.5)./s(1,t)).^2)+e;
%     aconslow = optimconstr(N-1,T);
%     aconslow(n,t)=a(n,t)>=alpha*(1-(v(n,t)./v0).^4-((s0+v(n,t).*hd+v(n,t).*(v(n,t)-v(n-1,t))/2/(alpha*beta).^0.5)./s(n,t)).^2)-e;
%     aconslow(1,t)=a(1,t)>=alpha*(1-(v(1,t)./v0).^4-((s0+v(1,t).*hd+v(1,t).*(v(1,t)-v(N,t))/2/(alpha*beta).^0.5)./s(1,t)).^2)-e;
idmcons = optimconstr(N-AVN,T);
idmcons(2:N-AVN,t)=a(HV,t)==alpha*(1-(v(HV,t)./v0).^4-((s0+v(HV,t).*hd+v(HV,t).*(v(HV,t)-v(HV-1,t))/2/(alpha*beta).^0.5)./s(HV,t)).^2);
idmcons(1,t)=a(1,t)==alpha*(1-(v(1,t)./v0).^4-((s0+v(1,t).*hd+v(1,t).*(v(1,t)-v(N,t))/2/(alpha*beta).^0.5)./s(1,t)).^2);

%dynamic for gap +(a(n-1,t-1)-a(n,t-1))*(ts.^2)/2
%+(a(N,t-1)-a(1,t-1))*(ts.^2)/2 +(init_a(n-1)-init_a(n))*(ts.^2)/2 +(init_a(N)-init_a(1))*(ts.^2)/2
gcons = optimconstr(N,T);
n=2:N;t=2:T;
gcons(n,t)=s(n,t)==s(n,t-1)+(v(n-1,t-1)-v(n,t-1))*ts+(a(n-1,t-1)-a(n,t-1))*(ts.^2)/2;
gcons(1,t)=s(1,t)==s(1,t-1)+(v(N,t-1)-v(1,t-1))*ts+(a(N,t-1)-a(1,t-1))*(ts.^2)/2;
gcons(n,1)=s(n,1)==init_s(n)+(init_v(n-1)-init_v(n))*ts+(init_a(n-1)-init_a(n))*(ts.^2)/2;
gcons(1,1)=s(1,1)==init_s(1)+(init_v(N)-init_v(1))*ts+(init_a(N)-init_a(1))*(ts.^2)/2;
lcons = optimconstr(T);
lcons(1:T) = sum(s(:,1:T))==Rlen-N*Vlen;
%dynamic for speed
vcons = optimconstr(N,T);nv=[1,HV];
vcons(nv,2:T)=v(nv,2:T)==v(nv,1:T-1)+a(nv,1:T-1)*ts+e1(nv,1:T-1)-e2(nv,1:T-1);
vcons(AV,2:T)=v(AV,2:T)==v(AV,1:T-1)+a(AV,1:T-1)*ts;
vcons(1:N,1)=v(1:N,1)==init_v(1:N)+init_a(1:N)*ts;
%dynamic for acceleration
jmax=2;
acons1 = optimconstr(AVN,T);
acons2 = optimconstr(AVN,T);
acons1(1:AVN,2:T)=a(AV,2:T)<=a(AV,1:T-1)+jmax*ts;
acons1(1:AVN,1)=a(AV,1)<=init_a(AV)+jmax*ts;
acons2(1:AVN,2:T)=a(AV,2:T)>=a(AV,1:T-1)-jmax*ts;
acons2(1:AVN,1)=a(AV,1)>=init_a(AV)-jmax*ts;

% % Monotonically increasing
% fcons = optimconstr(1);
% fcons(1) = sum(v,'all')>=result*1.1;

% ringroad.Constraints.fcons = fcons;
ringroad.Constraints.idmcons = idmcons;
ringroad.Constraints.gcons = gcons;
ringroad.Constraints.lcons = lcons;
ringroad.Constraints.vcons = vcons;
ringroad.Constraints.acons1 = acons1;
ringroad.Constraints.acons2 = acons2;
% -lamda*sum((a(n,t)-alpha*(1-(v(n,t)./v0).^4-((s0+v(n,t).*hd+v(n,t).*(v(n,t)-v(n-1,t))/2/(alpha*beta).^0.5)./s(n,t)).^2)).^2,'all').^0.5
%sum(v,'all')   -100*(sum(e1,"all")+sum(e2,'all'))
ringroad.Objective =sum(v,'all')-(sum(e1,"all")+sum(e2,'all'))*120/3.6;
end