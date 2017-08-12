cnt = 1;
for kkk = 123:134

clearvars -except SSS kkk cnt

%% Dataset selection

main_datasets_folder = '/home/sskim/Documents/bag_file';
dataset = 'sim';
dataset = strcat(dataset,num2str(kkk));
bag_file_name_bag = strcat(main_datasets_folder, '/', dataset);

%% Load the dataset
if exist(strcat(bag_file_name_bag, '.mat'))
    load(strcat(bag_file_name_bag, '.mat'))
else
    bag_file = strcat(bag_file_name_bag, '.bag');
    if exist(bag_file)
        bag_file_mat = strcat(main_datasets_folder, dataset, '.bag');
        data = convert_bag_to_mat(bag_file);
        save(strrep(bag_file_name_bag, '.bag', '.mat'), 'data')
    else
        disp('Error: dataset not found!');
        return;
    end
end

%%
time = data.hummingbird.flight_controller.feedback.time;

despos = data.hummingbird.flight_controller.feedback.desired_state.position;
pos = data.hummingbird.flight_controller.feedback.state_estimate.position;
desvel = data.hummingbird.flight_controller.feedback.desired_state.velocity;
vel = data.hummingbird.flight_controller.feedback.state_estimate.velocity;
desacc = data.hummingbird.flight_controller.feedback.desired_state.acceleration;
att = data.hummingbird.flight_controller.feedback.state_estimate.orientation;

%%
kp = diag([10 10 15]); kd = diag([4 4 6]);
ep = despos-pos;
ev = desvel-vel;
desThrust = kp*ep+kd*ev+desacc+[0;0;9.8];
e3 = [0 0 1]';
Se3 = [0 -1 0;1 0 0;0 0 0];

for k=1:size(desThrust,2)
   desDir(:,k) = desThrust(:,k)/norm(desThrust(:,k));
   q = att(4,k); q_ = att(1:3,k); q1 = q_(1); q2 = q_(2); q3 = q_(3);
   Sq_ = [ 0 -q3  q2;
          q3   0 -q1;
         -q2  q1   0];
   dir(:,k) = q*(q*e3-Se3*q_)+(e3'*q_)*q_+Sq_*(q*e3-Se3*q_);
   ang(k) = norm(cross(desDir(:,k),dir(:,k)));
end

%%
windowSize = 20;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
new_ang = filter(b,a,ang);

%%
for k=1:size(desvel,2)
   velnorm(k) = norm(desvel(:,k)); 
end
[~,idx1] = find(abs(velnorm(1,:)) < 1e-5);
[~,idx] = find(diff(idx1) > 2);
idx = idx1(idx)+10;

p1 = [2.2385 0.15;0.15 0.1737];

hold on
jj = 1;
for kk = 1:length(idx)-1
    ep = despos(jj,idx(kk))-pos(jj,idx(kk));
    ev = desvel(jj,idx(kk))-vel(jj,idx(kk));
    e = [ep;ev];
    
    maxAngErr = max(new_ang(idx(kk):idx(kk+1)-1));
%     if and(max(new_ang(idx(kk):idx(kk+1)-1)) < 0.5, and(e'*p1*e < 1.1, abs(ev) < 2))
%     if and(or(maxAngErr < 0.32, and(maxAngErr > 0.34, maxAngErr < 0.4)), and(1, abs(ev) < 2))
    if and(maxAngErr < 0.4, and(1, abs(ev) < 2))
        for k=1:200
            p = SSS{1}(:,k);
            p1 = [p(1) p(3);p(3) p(5)];
            if e'*p1*e < 1
                break;
            end
        end
        
        ts_pos = timeseries((despos(jj,idx(kk):idx(kk+1)-10)-pos(jj,idx(kk):idx(kk+1)-10))',time(idx(kk):idx(kk+1)-10));
        rs_pos = resample(ts_pos,time(idx(kk)):0.05:time(idx(kk+1)-10));
        ts_vel = timeseries((desvel(jj,idx(kk):idx(kk+1)-10)-vel(jj,idx(kk):idx(kk+1)-10))',time(idx(kk):idx(kk+1)-10));
        rs_vel = resample(ts_vel,time(idx(kk)):0.05:time(idx(kk+1)-10));
        
        figure(100)
        plot3(rs_pos.Time-rs_pos.Time(1)+0.05*(k+1),rs_pos.Data,rs_vel.Data,'linewidth',2);
    else
        cnt = cnt+1;
    end
end

end

view(-90,0)
% axis equal
axis([0 5 -1 1 -2.5 2.5])

