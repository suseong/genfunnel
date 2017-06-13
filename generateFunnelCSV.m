function generateFunnelCSV(time,cx,cy,cz,filename)

traj_matrix = [time cx cy cz];

csvwrite(filename,traj_matrix);

end