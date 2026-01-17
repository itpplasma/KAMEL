classdef managehdf5 < handle & hdf5_output
%##########################################################################
%classdef managehdf5 < handle & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to manage the reading and writing to hdf5 file.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path, fname
% READONLY:
% *)
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = managehdf5(path, fname)
% *) createfile(obj, overwrite)
% *) writetofile(obj, data) <- need to figure this out
% *) writegfiletohdf5(obj, path)
% *) writecfiletohdf5(obj, path)
% *) writeprofiletohdf5(obj, prof, path)
% *) readfromfile(obj, data)
% *) deletefile(obj)
%##########################################################################

%author:   Markus Markl
%created:  12.02.2020


properties
    
    fname       % name of the hdf5 file
    path        % path to the hdf5 file
    fullpath    % full path to file
    fullpathwodot % full path w/o hdf5 ending
    
end


methods
    function obj = managehdf5(path, fname)
        %##########################################
        % function obj = managehdft(path, fname)
        %##########################################
        % description:
        %------------------------------------------
        % Constructor of the hdf5 manage class.
        %##########################################
        %------------------------------------------
        % input:
        %------------------------------------------
        % path      ... path to the file
        % fname     ... name of the file, w/o .hdf5 ending
        %##########################################
        
        obj.path = path;
        obj.fname = fname;
        obj.fullpath = [obj.path, obj.fname, '.hdf5'];
        obj.fullpathwodot = [obj.path, obj.fname]; 
    
    end
    
    function change_path_fname(obj, path, fname)
        %##########################################
        % function change_path_fname(obj, path, fname)
        %##########################################
        % description:
        % Changes the objects path and file name
        %------------------------------------------
        % input:
        %------------------------------------------
        % path      ... path to the file
        % fname     ... name of the file, w/o .hdf5 ending
        %##########################################
        
        obj.path = path;
        obj.fname = fname;
    end
    
    function createfile(obj, overwrite)
        %##########################################
        % function createfile(obj, overwrite)
        %##########################################
        % description:
        % Method that creates the file with option
        % to overwrite existing file.
        %------------------------------------------
        % input:
        %------------------------------------------
        % overwrite     ... if unequal 0, an existing
        % file will be overwritten.
        %##########################################
        
        
        if(overwrite) % if overwrite is unequal 0, an
            % existing file will be deleted
            fulldir = [obj.path, obj.fname, '.hdf5'];
            if(exist(fulldir, 'file')) % if file exists, delete
                delete(fulldir);
            end
            
            % Create file
        else
            warning('File exists, no overwrite')
        end
        
    end
    
    function writegfiletohdf5(obj, pathtogfile)
        %##########################################
        % function writegfiletohdf5(obj, path)
        %##########################################
        % description:
        % Writes data of g-file to hdf5 file
        %------------------------------------------
        % input:
        %------------------------------------------
        % pathtogfile   ... path to the gfile
        %##########################################
        gfiledata = importdata(pathtogfile, ' ', 1);
        try
            h5create(obj.fullpath, '/input/gfile', size(gfiledata.data));
        catch
           warning('gfile already exists, will overwrite') 
        end
            
        h5write(obj.fullpath, '/input/gfile', gfiledata.data);
        % To do: add descriptive attribute
        
    end
    
    function writecfiletohdf5(obj, pathtocfile)
        %##########################################
        % function writecfiletohdf5(obj, path)
        %##########################################
        % description:
        % Writes data of a coil file to hdf5 file
        %------------------------------------------
        % input:
        %------------------------------------------
        % pathtocfile   ... path to the coil file
        %##########################################
        cfiledata = importdata(pathtocfile);
        try
            h5create(obj.fullpath, '/input/cfile', size(cfiledata));
        catch
            warning('cfile already exists, will overwrite')
        end
        
        h5write(obj.fullpath, '/input/cfile', cfiledata);
        % To do: add descriptive attribute
    end
    
    
    function writeprofiletohdf5(obj, prof, path)
        %##########################################
        % function writeprofiletohdf5(obj, prof, path)
        %##########################################
        % description:
        % Writes profile data to hdf5 file.
        %------------------------------------------
        % input:
        %------------------------------------------
        % path   ... path to the profile data
        % prof   ... name of profile that is written to hdf5 file
        %##########################################
        profiledata = importdata(path);
        try
            h5create(obj.fullpath, ['/input/profiles/', prof], size(profiledata));
        catch
            warning(['profile ', prof, ' already exists'])
        end
        
        h5write(obj.fullpath, ['/input/profiles/', prof], profiledata);
        % To do: add descriptive attribute
    end
    
    
    function writefluxdatatohdf5(obj, fluxdatapath)
        %##########################################
        % function writefluxdatatohdf5(obj, fluxdatapath)
        %##########################################
        % description:
        % Writes flux data to hdf5 file.
        %------------------------------------------
        % input:
        %------------------------------------------
        % fluxdatapath  ... path to flux data folder
        %##########################################
        btorrbig = importdata([fluxdatapath, 'btor_rbig.dat']);
        try
            h5create(obj.fullpath, '/input/fluxdata/btor', size(btorrbig(1)));
        catch
            warning('fluxdata/btor already exists');
        end
        
        h5write(obj.fullpath, '/input/fluxdata/btor', btorrbig(1));
        
        try
            h5create(obj.fullpath, '/input/fluxdata/bigr', size(btorrbig(2)));
        catch
            warning('fluxdata/bigr already exists');
        end
        
        h5write(obj.fullpath, '/input/fluxdata/bigr', btorrbig(2));
        
        equilrqpsi = importdata([fluxdatapath,'equil_r_q_psi.dat'], ' ', 3);

        try
            h5create(obj.fullpath, '/input/fluxdata/equil_r_q_psi', size(equilrqpsi.data));
        catch
            warning('fluxdata/equil_r_q_psi already exists');
        end
        
        h5write(obj.fullpath, '/input/fluxdata/equil_r_q_psi', equilrqpsi.data);
        
    end
    
    function writedata(obj, group, data)
        %##########################################
        % function writefluxdatatohdf5(obj, fluxdatapath)
        %##########################################
        % description:
        % Writes data to hdf5 file under given group.
        % (Used to write relevant GPEC data to hdf5 file. The 
        % relevant data is the current at the resonant surfaces
        % determined by the mode numbers.)
        %------------------------------------------
        % input:
        %------------------------------------------
        % group  ... specifies group in hdf5 file
        % data   ... data array to be written
        %##########################################
        try
            h5create(obj.fullpath, group, size(data));
        catch
            warning('fluxdata/equil_r_q_psi already exists');
        end
        
        h5write(obj.fullpath, group, data);
    end
    
    
    
    
    
    
    
    function readfromfile(obj, data)
        
    end
    
    function deletefile(obj)
        
    end
    
    
end


end
