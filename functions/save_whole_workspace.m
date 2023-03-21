function [obj_out,file_path] = save_whole_workspace(obj, zip_file_path, original_log_folder, original_logfile)
% adapted from save_workspace, Andrew Toader, The Rockefeller University
fprintf('------------- SAVE THE WHOLE WORKSPACE ----------\n\n');
            
            obj.compress_results();
            file_path = fullfile(obj.P.log_folder,  ['final_',strrep(get_date(), ' ', '_'), '.mat']);
            log_file = obj.P.log_file; 
            if exist('original_logfolder', 'var')
               obj.P.log_folder = original_log_folder;  
            end
            if exist('original_logfile', 'var')
                obj.P.log_file = original_logfile; 
            end
            evalin('caller', sprintf('save(''%s'', ''neuron'', ''save_*'', ''show_*'', ''use_parallel'', ''with_*'', ''-v7.3''); ', file_path));
            try
                fp = fopen(log_file, 'a');
                fprintf(fp, '\n--------%s--------\n[%s]\bSave the current workspace into file \n\t%s\n\n', get_date(), get_minute(), file_path);
                fprintf('The current workspace has been saved into file \n\t%s\n\n', file_path);
                fp.close();
            end
            
            if exist('zip_file_path', 'var') && ~isempty(zip_file_path)
               [zip_dir, zip_name, ~] = fileparts(get_fullname(zip_file_path)); 
               zip([zip_name, '.zip'], {file_path, log_file}, zip_dir); 
               fprintf('The results and the log files were compresed into file \n%s\n', zip_file_path); 
            end

            
            obj_out = obj;