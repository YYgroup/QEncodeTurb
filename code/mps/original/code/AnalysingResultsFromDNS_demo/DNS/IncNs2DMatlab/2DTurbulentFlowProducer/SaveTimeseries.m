savefile = strcat(strcat(CaseName,num2str(N)),'.mat');
save(savefile, 'Ekin', 'Diss', 'dt', 'N', 'Time', 'ElapsedTime');
