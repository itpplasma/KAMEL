module h5mod
    USE hdf5_tools
    integer(HID_T) :: h5_id, group_id_1, group_id_2, group_id_3, dataset_id
    logical :: h5_exists_log
    integer :: mode_n, mode_m
    character(len=1024) :: path2inp ! path to hdf5 file with input data
    character(len=1024) :: path2time ! path to hdf5 file from which time evolved profiles are read
    character(len=1024) :: path2out ! path to hdf5 file where output is written
    character(len=1024) :: h5_mode_groupname
end module h5mod