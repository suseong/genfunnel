    function sendTrajMsg(obj, xtraj)
      % Check trajectory
      if ~strcmp(class(xtraj), 'PPTrajectory')
        error('DemoController:   Trajectory is not of type PPTrajectory.');
      end
      
      traj_msg = rosmessage('state_info_msgs/TrajectoryMatrix');
      
      % Get properties
      N = xtraj.pp.pieces;
      t_breaks = xtraj.pp.breaks;
      coefs = xtraj.pp.coefs;
      states = floor(size(coefs,1)/N);
      order = xtraj.pp.order;
      
      % Reconstruct matrix in a nice way
      traj_matrix = zeros(N, 1+order*states);
      for i=1:N
        entry = zeros(1,order*states);
        for j=1:states
          entry((1+order*(j-1)):(order*(j))) = coefs((i-1)*states+j,:);
        end
        traj_matrix(i,:) = [t_breaks(i+1) entry];
      end
      
      traj_serial = reshape(traj_matrix',[N*(1+order*states),1]);
      
      % Fill in message.
      traj_msg.Id.Data = 0;
      traj_msg.States.Data = states;
      traj_msg.Nodes.Data = N;
      traj_msg.Order.Data = order;
      [traj_msg.Rows.Data, traj_msg.Columns.Data] = size(traj_matrix);
      traj_msg.Time.Data = traj_matrix(end,1);
      traj_msg.Sequence.Data = traj_serial;
      
      % Send the message.
      send(obj.traj_pub, traj_msg);
    end
