function add_rm_paths(add_rm_flag)
% from https://github.com/daveckman/plausible-screening/blob/master/scripts/add_rm_paths.m

folders = strcat({'../'}, {'src','data','figures'});

switch add_rm_flag
    case 'add'
        for k = 1:length(folders)
            addpath(folders{k})
        end

    case 'remove'
        for k = 1:length(folders)
            rmpath(folders{k})
        end

    otherwise
        disp('\nERROR: Invalid argument for add_rm_paths(). Argument must be "add" or "remove".\n')
end

end

