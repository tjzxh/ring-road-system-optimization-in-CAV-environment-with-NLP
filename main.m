%% Warm up
% Create platoon array
total_time = 200;step = 1;N = 22;T=100;ts=1;Rlen=230;Vlen=5;
alpha=0.73;beta=1.67;v0=30;hd=1.5;s0=2;
spacing = zeros(total_time, N);speed=zeros(total_time, N); acce=zeros(total_time, N);
% Simulation start +2*rand()-1
for veh_num=1:N
    spacing(1,veh_num) = (Rlen-N*Vlen)/N;
    speed(1,veh_num) = 30/3.6;
    acce(1,veh_num) = 0;
end
% Update speed and position with kinematic equation
for time_step = 2:total_time
    for veh_num=1:N
        speed(time_step,veh_num) = speed(time_step-1,veh_num) + acce(time_step-1,veh_num) * step;
        if speed(time_step,veh_num) < 0
            speed(time_step,veh_num) = 0;
            acce(time_step-1,veh_num) = -speed(time_step-1,veh_num);
        end
        if veh_num==1
            % Ring road
            % +(acce(time_step-1,N)-acce(time_step-1,veh_num))*(step.^2)/2 +(acce(time_step-1,veh_num-1)-acce(time_step-1,veh_num))*(step.^2)/2
            spacing(time_step,veh_num) = spacing(time_step-1,veh_num) + (speed(time_step-1,N)-speed(time_step-1,veh_num)) * step+(acce(time_step-1,N)-acce(time_step-1,veh_num))*(step.^2)/2;
            vd = speed(time_step,veh_num) - speed(time_step,N);
        else
            spacing(time_step,veh_num) = spacing(time_step-1,veh_num) + (speed(time_step-1,veh_num-1)-speed(time_step-1,veh_num)) * step+(acce(time_step-1,veh_num-1)-acce(time_step-1,veh_num))*(step.^2)/2;
            vd = speed(time_step,veh_num) - speed(time_step,veh_num-1);
        end
        % Calculate acceleration with IDM
        acce(time_step,veh_num)= alpha*(1-(speed(time_step,veh_num)/v0).^4-((s0+speed(time_step,veh_num)*hd+speed(time_step,veh_num)*vd/2/(alpha*beta).^0.5)/spacing(time_step,veh_num)).^2);
        if acce(time_step,veh_num)<-2
            acce(time_step,veh_num)=-2;
        end
    end
end
% First optimization
init_a=acce(total_time-T-1,:)';
init_v=speed(total_time-T,:)';
init_s=spacing(total_time-T,:)';
ringroad = rrso(init_a,init_v,init_s,0);
init.a =acce(total_time-T:total_time-1,:)';init.v =speed(total_time-T+1:total_time,:)';init.s =spacing(total_time-T+1:total_time,:)';
init.e1=zeros(N,T);init.e2=zeros(N,T);
opts = optimoptions('fmincon','Display','iter',"EnableFeasibilityMode",true,"SubproblemAlgorithm","cg");
% opts = optimoptions('linprog','Display','iter',"ConstraintTolerance",0.001);
[sol,fval,exitflag,output] = solve(ringroad,init,"Options",opts);
% % Loop for optimization
% iter_time = 1;
% results = zeros(10,1);
% while iter_time <= 10
%     init_a=sol.a(:,T-1);
%     init_v=sol.v(:,T);
%     init_s=sol.s(:,T);
%     init.a=sol.a;init.v=sol.v;init.s=sol.s;init.e1=zeros(N,T);init.e2=zeros(N,T);
%     ringroad = rrso(init_a,init_v,init_s,sum(sol.v,"all"));
%     [sol,fval,exitflag,output] = solve(ringroad,init,'Options',opts);
%     results(iter_time)=sum(sol.v,"all")/N/T;
%     iter_time = iter_time + 1;
% end