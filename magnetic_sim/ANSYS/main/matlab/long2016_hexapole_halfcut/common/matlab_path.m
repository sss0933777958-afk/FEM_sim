function p = matlab_path(model, varargin)
% DEPRECATED (2026-06-26): MATLAB_data/ 已全量遷移到各活動的 local data/ 並移除
%   （規則 rules/matlab-output-layout.md）。`.mat` 一律放 matlab/<model>/<activity>/data/，
%   不再經本 resolver。本檔保留僅供歷史/相容；其回傳路徑已不存在，請勿再用於新程式。
%MATLAB_PATH  Central resolver for the design (magnetic_sim/ANSYS/main) MATLAB-output root (.mat/.csv/.npz).
%   Resolved RELATIVE to this file (no hard-coded drive letter). MATLAB analysis
%   outputs are organised as magnetic_sim/ANSYS/main/MATLAB_data/<model>/<function>/ (function =
%   charge_fit / bs_matrix / flux_profile / freq_response / bh_saturation / ...).
%
%   p = matlab_path()                          -> .../main/MATLAB_data
%   p = matlab_path(model)                     -> .../main/MATLAB_data/<model>
%   p = matlab_path(model, 'charge_fit')       -> .../<model>/charge_fit
%   p = matlab_path(model, 'charge_fit', f)    -> .../<model>/charge_fit/<f>
%
%   FEM data (.dat/.db) lives under magnetic_sim/ANSYS/main/ANSYS_data/ -> use ansys_path() for those.
%   If this folder is ever renamed, change only MAT_DIRNAME below.

    MAT_DIRNAME = 'MATLAB_data';                        % <== the ONLY place the name lives

    here = fileparts(mfilename('fullpath'));            % .../main/matlab/long2016_hexapole_halfcut/common
    root = fileparts(fileparts(fileparts(here)));       % .../main
    base = fullfile(root, MAT_DIRNAME);                 % .../main/MATLAB_data

    if nargin == 0
        p = base;
    else
        p = fullfile(base, model, varargin{:});
    end
end
