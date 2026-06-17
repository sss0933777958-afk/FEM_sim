function p = ansys_path(model, varargin)
%ANSYS_PATH  Central resolver for the design (magnetic_sim/ANSYS/main) FEM/fitting data root.
%   Returns an absolute path under the main data tree, resolved RELATIVE to this
%   file's own location (no hard-coded drive letter), so the repo can be moved
%   and the data-folder name can be renamed by editing ONE line here.
%
%   p = ansys_path()                      -> .../main/ANSYS_data
%   p = ansys_path(model)                 -> .../main/ANSYS_data/<model>
%   p = ansys_path(model, 'coil1','standard') -> .../main/ANSYS_data/<model>/coil1/standard
%
%   FEM data (.dat/.db) lives under magnetic_sim/ANSYS/main/ANSYS_data/. MATLAB outputs (.mat) live under
%   magnetic_sim/ANSYS/main/MATLAB_data/<model>/<function>/  -> use matlab_path() for those.
%
%   Usage in scripts (replaces hard-coded 'G:\...\main\ANSYS_data\<model>'):
%       results_root = ansys_path('long2016_hexapole_halfcut');
%       d = import_ansys_data(ansys_path(model,'coil1','standard'), 'wp', 'coil1');
%
%   If the FEM data folder is ever renamed again, change only DATA_DIRNAME below.

    DATA_DIRNAME = 'ANSYS_data';                        % <== the ONLY place the name lives

    here = fileparts(mfilename('fullpath'));            % .../main/matlab/long2016_hexapole_halfcut/common
    root = fileparts(fileparts(fileparts(here)));       % .../main
    base = fullfile(root, DATA_DIRNAME);                % .../main/ANSYS_data

    if nargin == 0
        p = base;
    else
        p = fullfile(base, model, varargin{:});
    end
end
