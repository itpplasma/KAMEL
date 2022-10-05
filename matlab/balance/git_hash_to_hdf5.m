function git_hash_to_hdf5(h5out)

    s = which(mfilename);
    s = [fileparts(s), '/../'];
    s = what(s);
    path = [s.path, '/'];

    [status, git_hash_string] = system(['git -C ', path, ' rev-parse HEAD']);
    dim = size(git_hash_string, 1);


    
    file_id = H5F.open(h5out, 'H5F_ACC_RDWR', 'H5P_DEFAULT')
    filetype = H5T.copy('H5T_C_S1')
    H5T.set_size(filetype, dim);
    memtype = H5T.copy('H5T_C_S1');
    H5T.set_size(memtype, dim);

    space_id = H5S.create_simple(1, dim, []);
    dataset_id = H5D.create(file_id, 'String', filetype, space_id, 'H5P_DEFAULT');

    H5D.write(dataset_id, memtype, 'H5S_ALL', 'H5P_DEFAULT', git_hash_string);

    H5D.close(dataset_id);
    H5S.close(space_id);
    H5T.close(filetype);
    H5T.close(memtype);
    H5F.close(file_id);

%    h5write(h5out, '/git', git_hash_string)

end