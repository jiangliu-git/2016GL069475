themis_folder = '/Weiyun/data/themis'
ground_folder = '/Weiyun/data/ground'
sw_folder = '/Weiyun/data/sw'
plasma_folder = '/plasma_fill'

if keyword_set(computer) then begin
	fgs_folder=computer + themis_folder + '/fgs'
	fgl_folder=computer + themis_folder + '/fgl'
	pos_folder=computer + themis_folder + '/pos'
	efs_folder=computer + themis_folder + '/efs_dsl'
	eff_folder=computer + themis_folder + '/eff_dot0_dsl'
	Pall_folder=computer + themis_folder + plasma_folder + '/Pall'
	beta_folder=computer + themis_folder + plasma_folder + '/beta'
	Pth_folder=computer + themis_folder + plasma_folder + '/Pth'
	vi_folder = computer + themis_folder + plasma_folder + '/vi'
	viperp_folder = computer + themis_folder + plasma_folder + '/viperp'
	ve_folder = computer + themis_folder + plasma_folder + '/ve'
	vexb_folder = computer + themis_folder + plasma_folder + '/vexb'
	veperp_folder = computer + themis_folder + plasma_folder + '/veperp'
	Pttl_fgs_folder=computer+themis_folder+plasma_folder+'/Pttl_fgs'
	Blobe_folder=computer+themis_folder+plasma_folder+'/Blobe'
	al_folder = computer+ground_folder+'/pseudo_al'
	kyoto_al_folder = computer+ground_folder+'/kyoto_al'
	omni_b_folder = computer+sw_folder+'/omni_b_gsm'
endif

list_folder='../lists'
save_folder='variables'
pic_folder = '../../../results/temp'

;;;;;; 
l_thick = 2.4
;;; constants and signs
Re = 6371.
mu0 = 4*!pi*1e-7
nTesla2_to_nPa = 0.01/25.132741
perp_sign = '!9'+string("136B)+'!X'
cross = '!9'+string("264B)+'!X'
alpha_letter = '!9'+string("141B)+'!X'
beta_letter = '!9'+string("142B)+'!X'
theta_letter = '!9'+string("161B)+'!X'
phi_letter = '!9'+string("152B)+'!X'
cap_theta_letter = '!9'+string("121B)+'!X'
delta_letter = '!9'+string("144B)+'!X'
cap_delta_letter = '!9'+string("104B)+'!X'
jiao_l = '!9'+string("341B)+'!X'
jiao_r = '!9'+string("361B)+'!X'
dot_sign = '!9'+string("327B)+'!X'
plus_sign = '!9'+string("53B)+'!X'
minus_sign = '!9'+string("55B)+'!X'
angle_sign = '!9'+string("320B)+'!X'

;;; for DFB properties
seconds_check = 15.

;;; turn off the timestamp of tplot
time_stamp,/off

;;; for panels layout
left_margin = 0.2
right_margin = 0.07
top_margin = 0.05
bot_margin = 0.07
space_horiz = 0.1
space_vert = 0.01
