function write_version_to_hdf5(h5out, version_number)

    s = which(mfilename);
    s = [fileparts(s), '/../'];
    s = what(s);
    path = [s.path, '/'];

    [status, git_hash_string] = system(['git -C ', path, ' rev-parse HEAD']);
    dim = size(git_hash_string, 1);


    disp(' ')
    disp(['Git commit hash: ', git_hash_string])
    disp(' ')
    h5create(h5out, '/code_version', 1);
    h5write(h5out, '/code_version', version_number);
    h5writeatt(h5out, '/code_version', 'git_hash', git_hash_string);

end