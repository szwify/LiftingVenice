function [outputArg1,outputArg2] = save_to_csv(u_tot,name)
for t=1:n_step
    dir_path = '\paraview';
    file_name = sprintf('test.csv.%d',[t]);
    out = fullfile(dir_path,file_name);
    fileID = fopen(out,'w');
    fprintf(fileID,'X, Y, norm \n');
    fprintf(fileID,'%d, %d, %d \n',reshape([coordonnees,norm_u]',[],1)');
    fclose(fileID);
end
end

