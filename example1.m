clear all

%%----------Download---------
% Magnetosheath
%tint = [irf_time([2002 03 27 10 16 00]) irf_time([2002 03 27 10 17 30])];
% Magnetopause
tint = [irf_time([2002 02 20 13 21 00]) irf_time([2002 02 20 13 23 00])];
% Magnetotail
%tint = [irf_time([2003 08 24 18 43 00]) irf_time([2003 08 24 18 44 30])];

%B-field
caa_download(tint,'C?_CP_FGM_FULL');


gseMagC1 = irf_get_data('B_vec_xyz_gse__C1_CP_FGM_FULL','caa','mat');
gseMagC2 = irf_get_data('B_vec_xyz_gse__C2_CP_FGM_FULL','caa','mat');
gseMagC3 = irf_get_data('B_vec_xyz_gse__C3_CP_FGM_FULL','caa','mat');
gseMagC4 = irf_get_data('B_vec_xyz_gse__C4_CP_FGM_FULL','caa','mat');

%Downloads position data.
caa_download(tint,'C?_CP_AUX_POSGSE_1M')
R = [];
R.R1 = irf_get_data('sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M','caa','mat');
R.R2 = irf_get_data('sc_r_xyz_gse__C2_CP_AUX_POSGSE_1M','caa','mat');
R.R3 = irf_get_data('sc_r_xyz_gse__C3_CP_AUX_POSGSE_1M','caa','mat');
R.R4 = irf_get_data('sc_r_xyz_gse__C4_CP_AUX_POSGSE_1M','caa','mat');



%calls c_4_timing_mva for the z-axis
c_4_v_timing_mva('gseMagC?',R,4);