function genCSV(time,cx,cy,cz,filename)

traj_matrix = [time; cx; cy; cz];
dlmwrite(filename,traj_matrix,'delimiter', ',' , 'precision', 15);

end