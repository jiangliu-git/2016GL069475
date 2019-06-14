pro main_current_stat
;; generate bin plots, distribution of ctcs, FACs in the tail.
;; load list, instead of seasons as in main_ctcs_stat.
thm_init

computer='I:'

@folders

;;; choose list
list_suf = '' ;; all times
;list_suf = '_grw_csv' ;; conservative growth phase
;list_suf = '_sub' ;; substorm
;list_suf = '_exp_csv' ;; conservative (early) expansion phase
;list_suf = '_rcv_lbl' ;; liberal recovery phase
;list_suf = '_grwexp' ;; late growth (5min) + early expansion phase
;list_suf = '_grwexp20' ;; late growth (20min) + early expansion phase
;list_suf = '_nongrwexp' ;; every but grwexp

;;; choose method
;method_suf = '' ;; use seperated probes in Z
method_suf = '_3p' ;; use three probes to compute j in the direction perp to the plane formed by the three

;;;;;; load lists
if strcmp(list_suf, '') then begin
	;; all
	season1 = ['2007 11 1','2008 6 30']
	season2 = ['2008 12 1','2009 6 30']
	season3 = ['2010 2 1','2010 6 30']
	season4 = ['2011 2 1','2011 7 31']
	season5 = ['2012 4 13','2012 10 14']
	season6 = ['2013 6 15','2013 10 15']
	season7 = ['2014 6 29','2014 10 16']
	;;;; seasons ;;
	n_split = 40
	seasons = [[season1],[season2], [season3], [season4], [season5], [season6], [season7]]
	;;; split the seasons to a list, each split to n_split time ranges
	events = strarr(2, n_split*n_elements(seasons[0,*]))
	for i = 0, n_elements(seasons[0,*])-1 do begin
		season = time_double(seasons[*,i])
		times = linspace(season[0], season[1], n_split+1)
		tranges_season = [transpose(times[0:n_split-1]), transpose(times[1:n_split])]
		events[*, i*n_split:(i+1)*n_split-1] = time_string(tranges_season)
	endfor
endif else begin
	;; lists
	case list_suf of
	'_sub': listname = 'sb_substorm_list'
	else: listname = 'sb'+list_suf+'_list'
	endcase
	
	events = load_trange_list(list_folder+'/'+listname+'.txt')
endelse

probes_more = ['a','b','c','d','e']
probes_less = ['a','c','d']

;;;;; settings for fast processing
reduce_pts = 1 ;; reduce the position data points to this ratio, since pos is 1-min res, this represents minutes.
events_sep = 10. ;; in minutes, to tell events apart

;;; range to search
xrange = [-30.,-6.]
rhomax = 12.

;;; sc requirements to compute j

;; 2p
dsep_max = 2. ;; in RE, need smaller than this
dother_ratio2sep = 1. ;; in ratio, need smaller than this. Must be smaller than 1, or every thing will be wrong

;; 3p
dsep_max_3p = 2. ;; in RE, need smaller than this
angle_3p_pos = 15. ;; in degrees, to tell whether the 3-p probes are in a good triangle. The angles of the triangle should not be smaller than this.
;angle_3p = 30. ;; in degrees, to tell whether the 3-p computed j is close enough to a direction.
;angle_3p = 45. ;; in degrees, to tell whether the 3-p computed j is close enough to a direction.
angle_3p = 90. ;; disable angle_3p by setting it to 90


;;; result arrays for increment
i_events = 0
t_multi = 0d ;; time
pos_ave = transpose([0., 0., 0.])
beta_ave = 0.
beta_far = 0.
beta_near = 0.
h_ave = 0.
h_far = 0.
h_near = 0.
hxy_ave = 0.
hxy_far = 0.
hxy_near = 0.
al_compare = 0. ;; this is pseudo AL
al_kyoto = 0. ;; this is kyoto AL
l_harris = 0.
z_sigma = 0.
l_sigma = 0.
zns = 0.
zns_far = 0.
zns_near = 0.
Blobe = 0.
b_ave = transpose([0., 0., 0.])
imf = transpose([0., 0., 0.])
if ~strcmp(method_suf, '_3p') then begin
	;;;; two probes only
	jy = 0.
	jx = 0.
	dxyz = transpose([0., 0., 0.])
	probes_dual = transpose(['', ''])
endif else begin
	;;; three probes only
	dir_3p_type = 0. ;; 1: close to X; 2: close to Bdir; 3: close to both
	angle_3p_jvx = 0. ;; in degs
	angle_3p_jvb = 0. ;; in degs
	j3p_regx = 0. ;; current value regarding X direction
	j3p_regb = 0. ;; current value regarding B direction
	j3p = transpose([0., 0., 0.])
	distance_3p = transpose([0., 0., 0.])
	triangle_3p = transpose([0., 0., 0.])
	probes_3p = transpose(['', '', ''])
endelse

;;; mark the probes
;probes_dual = 'probes'

for i = 0, n_elements(events[0,*])-1 do begin
	trange_load = time_double(events[*,i])
	del_data, 'th?_state_pos*'
	load_bin_data, trange = trange_load, probes = probes, datatype = 'pos', datafolder = pos_folder, /tclip
	del_data, 'th?_fgs_gsm*'
	load_bin_data, trange = trange_load, probes = probes, datatype = 'fgs', datafolder = fgs_folder, /tclip
	del_data, 'th?_beta*'
	load_bin_data, trange = trange_load, probes = probes, datatype = 'beta', datafolder = beta_folder, /tclip
	del_data, 'thg_idx_al*'
	load_bin_data, trange = trange_load, probes = probes, datatype = 'pseudo_al', datafolder = al_folder, /tclip
	del_data, 'kyoto_al*'
	load_bin_data, trange = trange_load, probes = probes, datatype = 'kyoto_al', datafolder = kyoto_al_folder, /tclip
	del_data, 'th?_Pall*'
	load_bin_data, trange = trange_load, probes = probes, datatype = 'Pall', datafolder = Pall_folder, /tclip
	del_data, 'omni_b_gsm*'
	load_bin_data, trange = trange_load, datatype = 'omni_b_gsm', datafolder = omni_b_folder, /tclip
	;; load mode to tell it's in fast survey
	;del_data, 'th?_hsk_issr_mode_raw*'
	;load_bin_data, trange = trange_load, probes = probes, datatype = 'mode', datafolder = mode_folder, /tclip
	;;;; generate the time arrays and interpolate all quantities to this
	n_t_arr = ceil((trange_load[1]-trange_load[0])/(60.*reduce_pts))
	if n_t_arr eq 0 then n_t_arr = 1
	t_all = linspace(trange_load[0], trange_load[1], n_t_arr)
	pos_all = fltarr(n_t_arr, n_elements(probes), 3)+!values.f_nan
	beta_all = fltarr(n_t_arr, n_elements(probes))+!values.f_nan
	Pth_all = fltarr(n_t_arr, n_elements(probes))+!values.f_nan
	b_all = fltarr(n_t_arr, n_elements(probes), 3)+!values.f_nan
	;; true or false arrays
	pos_inrange = bytarr(n_t_arr, n_elements(probes))
	;;;;;;; interpolate
	;; single quantities
	if tv_exist('thg_idx_al_tclip') then begin
		get_data, 'thg_idx_al_tclip', t_al, al
		if n_elements(t_al) gt 1 then begin
			al_all = interpol(al, t_al, t_all)
			;; mark out of range points to be NaN or strange results will appear.
			i_out = where((t_all lt t_al[0]-reduce_pts*60.) or (t_all gt t_al[-1]+reduce_pts*60.), n_out)
			if n_out gt 0 then al_all[i_out] = !values.f_nan
		endif else al_all = replicate(al, n_t_arr)
	endif else begin
		al_all = fltarr(n_t_arr)+!values.f_nan ;; al index
	endelse
	if tv_exist('kyoto_al_tclip') then begin
		get_data, 'kyoto_al_tclip', t_kal, kal
		if n_elements(t_kal) gt 1 then begin
			kal_all = interpol(kal, t_kal, t_all)
			;; mark out of range points to be NaN or strange results will appear.
			i_out = where((t_all lt t_kal[0]-reduce_pts*60.) or (t_all gt t_kal[-1]+reduce_pts*60.), n_out)
			if n_out gt 0 then kal_all[i_out] = !values.f_nan
		endif else kal_all = replicate(kal, n_t_arr)
	endif else begin
		kal_all = fltarr(n_t_arr)+!values.f_nan ;; al index
	endelse
	if tv_exist('omni_b_gsm_tclip') then begin
		get_data, 'omni_b_gsm_tclip', t_imf, b_imf
		if n_elements(t_imf) gt 1 then begin
			imf_all = fltarr(n_t_arr, 3)
			for i_comp = 0, 2 do begin
				imf_all[*,i_comp] = interpol(b_imf[*,i_comp], t_imf, t_all)
			endfor
			;; mark out of range points to be NaN or strange results will appear.
			i_out = where((t_all lt t_imf[0]-reduce_pts*60.) or (t_all gt t_imf[-1]+reduce_pts*60.), n_out)
			if n_out gt 0 then imf_all[i_out,*] = !values.f_nan
		endif else imf_all = rebin(b_imf, n_t_arr, 3)
	endif else begin
	stop
		imf_all = fltarr(n_t_arr, 3)+!vimfues.f_nan ;; imf index
	endelse
	;; probe-dependent quantities
	for k = 0, n_elements(probes)-1 do begin
		sc = probes[k]
		if tv_exist('th'+sc+'_state_pos_tclip') then begin
		  	get_data, 'th'+sc+'_state_pos_tclip', t_pos, xyz
			xyz_re = xyz/RE
			for i_comp = 0, 2 do begin
				if n_elements(t_pos) gt 1 then begin
					pos_all[*,k,i_comp] = interpol(xyz_re[*,i_comp], t_pos, t_all)
					;; mark out of range points to be NaN or strange results will appear.
					i_out = where((t_all lt t_pos[0]-reduce_pts*60.) or (t_all gt t_pos[-1]+reduce_pts*60.), n_out)
					if n_out gt 0 then pos_all[i_out,k,i_comp] = !values.f_nan
				endif else pos_all[*,k,i_comp] = replicate(xyz_re[*,i_comp], n_t_arr) 
			endfor
			pos_inrange[*,k] = (pos_all[*,k,0] gt xrange[0]) and (pos_all[*,k,0] lt xrange[1]) and (sqrt(pos_all[*,k,1]^2+pos_all[*,k,2]^2) lt rhomax)
		endif
		if tv_exist('th'+sc+'_fgs_gsm_tclip') then begin
		  	get_data, 'th'+sc+'_fgs_gsm_tclip', t_b, bxyz
			for i_comp = 0, 2 do begin
				if n_elements(t_b) gt 1 then begin
					b_all[*,k,i_comp] = interpol(bxyz[*,i_comp], t_b, t_all)
					;; mark out of range points to be NaN or strange results will appear.
					i_out = where((t_all lt t_b[0]-reduce_pts*60.) or (t_all gt t_b[-1]+reduce_pts*60.), n_out)
					if n_out gt 0 then b_all[i_out,k,i_comp] = !values.f_nan
				endif else b_all[*,k,i_comp] = replicate(bxyz[*,i_comp], n_t_arr)
			endfor
		endif
		if tv_exist('th'+sc+'_beta_tclip') then begin
		  	get_data, 'th'+sc+'_beta_tclip', t_beta, beta_this
			if n_elements(t_beta) gt 1 then begin
				beta_all[*,k] = interpol(beta_this, t_beta, t_all)
				;; mark out of range points to be NaN or strange results will appear.
				i_out = where((t_all lt t_beta[0]-reduce_pts*60.) or (t_all gt t_beta[-1]+reduce_pts*60.), n_out)
				if n_out gt 0 then beta_all[i_out,k] = !values.f_nan
			endif else beta_all[*,k] = replicate(beta_this, n_t_arr)
		endif
		if tv_exist('th'+sc+'_Pall_tclip') then begin
		  	get_data, 'th'+sc+'_Pall_tclip', t_p, Pall
			if n_elements(t_p) gt 1 then begin
				Pth_all[*,k] = interpol(Pall[*,1], t_p, t_all)
				;; mark out of range points to be NaN or strange results will appear.
				i_out = where((t_all lt t_p[0]-reduce_pts*60.) or (t_all gt t_p[-1]+reduce_pts*60.), n_out)
				if n_out gt 0 then Pth_all[i_out,k] = !values.f_nan
			endif else Pth_all[*,k] = replicate(Pall[*,1], n_t_arr)
		endif
	endfor ;; for of k, probes

	;;; compute h
	Pttl_2 = nTesla2_to_nPa*(b_all[*,*,0]^2+b_all[*,*,1]^2)+Pth_all
	Bl = sqrt(Pttl_2/nTesla2_to_nPa)
	h_all = b_all[*,*,0]/Bl ;; scale height
	hxy_all = sqrt(b_all[*,*,0]^2+b_all[*,*,1]^2)/Bl ;; scale height

	;;; find suitable ranges and save the current
	for i_sc = 0, n_elements(probes)-2 do begin
		for j_sc = i_sc+1, n_elements(probes)-1 do begin
			d_pos = reform(pos_all[*,i_sc,*]-pos_all[*,j_sc,*], n_t_arr, 3)
			if strcmp(method_suf, '_3p') and n_elements(probes) ge 3 then begin
				;;;;;;; compute from 3 probes
				for k_sc = j_sc+1, n_elements(probes)-1 do begin
					;;; compute the current vector no matter what, this is a vector, in the direction of the plane formed by the three probes (make it 3xn_t_arr)
					j_this = curlometer3p(transpose(reform(b_all[*,i_sc,*])), transpose(reform(b_all[*,j_sc,*])), transpose(reform(b_all[*,k_sc,*])), $
								transpose(reform(pos_all[*,i_sc,*])), transpose(reform(pos_all[*,j_sc,*])), transpose(reform(pos_all[*,k_sc,*])), /RE, $
								distance_12 = distance_ij, distance_23 = distance_jk, distance_31 = distance_ki, position_mean = pos_mean, $
								angle_12 = angle_ij, angle_23 = angle_jk, angle_31 = angle_ki)

					;;; decide whether the satellite positions are good and 
					tell_good_pos = pos_inrange[*,i_sc] and pos_inrange[*,j_sc] and pos_inrange[*,k_sc];; good location
					tell_good_sep = (distance_ij lt dsep_max_3p) and (distance_jk lt dsep_max_3p) and (distance_ki lt dsep_max_3p)
					tell_good_triangle = (angle_ij gt angle_3p_pos) and (angle_jk gt angle_3p_pos) and (angle_ki gt angle_3p_pos)
					tell_good_config = tell_good_pos and tell_good_sep and tell_good_triangle

					;;; decide whether the computed current is closer to X or FAC locations
					b_ave_all = transpose(mean(b_all[*,[i_sc,j_sc,k_sc],*], dim=2, /nan)) ;; the average direction of the field lines of the three spacecraft
					dir_type_this = intarr(n_t_arr)
					angle_j_x = angle_vectors(j_this, rebin([1.,0,0], 3, n_t_arr), /deg)
					angle_j_b = angle_vectors(j_this, b_ave_all, /deg)
					;;; then store it if satisfy all.
					i_good_x = where(((angle_j_x lt angle_3p) or (angle_j_x gt 180.-angle_3p)) and tell_good_config, n_good_x)
					if n_good_x gt 0 then dir_type_this[i_good_x] = dir_type_this[i_good_x]+1
					i_good_b = where(((angle_j_b lt angle_3p) or (angle_j_b gt 180.-angle_3p)) and tell_good_config, n_good_b)
					if n_good_b gt 0 then dir_type_this[i_good_b] = dir_type_this[i_good_b]+2
					i_good = where(dir_type_this gt 0, n_good)
					if n_good gt 0 then begin
						;; compute the strength and sign of the current
						dir_type_this = dir_type_this[i_good]
						j_this = j_this[*, i_good]
						j_str = sqrt(total(j_this^2, 1))
						b_ave_this = b_ave_all[*, i_good]
						sign_j_x = float(transpose(sign(dotp_long(j_this, rebin([1.,0,0], 3, n_good)))))
						sign_j_b = float(transpose(sign(dotp_long(j_this, b_ave_this))*sign(b_ave_this[0,*])))
						i_bad_x = where((dir_type_this ne 1) and (dir_type_this ne 3), n_bad_x)
						if n_bad_x gt 0 then sign_j_x[i_bad_x] = !values.f_nan
						i_bad_b = where((dir_type_this ne 2) and (dir_type_this ne 3), n_bad_b)
						if n_bad_b gt 0 then sign_j_b[i_bad_b] = !values.f_nan

						;;;; add to all result arrays
						;; store two current values regarding two different directions
						j3p_regx = [j3p_regx, j_str*sign_j_x]
						j3p_regb = [j3p_regb, j_str*sign_j_b]
						j3p = [j3p, transpose(j_this)]

						;; store other data
						i_events = [i_events, replicate(i, n_good)]
						dir_3p_type = [dir_3p_type, dir_type_this]
						angle_3p_jvx = [angle_3p_jvx, angle_j_x[i_good]]
						angle_3p_jvb = [angle_3p_jvb, angle_j_b[i_good]]
						t_multi = [t_multi, t_all[i_good]]
						distance_3p_this = [[distance_ij[i_good]], [distance_jk[i_good]], [distance_ki[i_good]]]
						distance_3p = [distance_3p, distance_3p_this]
						triangle_3p_this = [[angle_ij[i_good]], [angle_jk[i_good]], [angle_ki[i_good]]]
						triangle_3p = [triangle_3p, triangle_3p_this]
						pos_these = pos_all[i_good,[i_sc,j_sc,k_sc],*]
						pos_ave = [pos_ave, mean(pos_these, dim=2, /nan)]
						b_ave = [b_ave, transpose(b_ave_this)]
						h_these = h_all[i_good,[i_sc,j_sc,k_sc],*]
						h_ave = [h_ave, mean(h_these, dim=2, /nan)]
						h_max_abs = max(abs(h_these), i_max_abs, dim=2, /nan)
						h_far = [h_far, h_max_abs*sign(h_these[i_max_abs])]
						h_min_abs = min(abs(h_these), i_min_abs, dim=2, /nan) ;; this is possibly smaller than h_ave because of signs when computing h_ave.
						h_near = [h_near, h_min_abs*sign(h_these[i_min_abs])]
						hxy_these = hxy_all[i_good,[i_sc,j_sc,k_sc],*]
						hxy_ave = [hxy_ave, mean(hxy_these, dim=2, /nan)]
						hxy_max_abs = max(abs(hxy_these), i_max_abs, dim=2, /nan)
						hxy_far = [hxy_far, hxy_max_abs*sign(hxy_these[i_max_abs])]
						hxy_min_abs = min(abs(hxy_these), i_min_abs, dim=2, /nan)
						hxy_near = [hxy_near, hxy_min_abs*sign(hxy_these[i_min_abs])]
						beta_ave = [beta_ave, abs(mean(beta_all[i_good,[i_sc,j_sc,k_sc],*]*sign(h_these), dim=2, /nan))]
						beta_far = [beta_far, min(beta_all[i_good,[i_sc,j_sc,k_sc],*], dim=2, /nan)]
						beta_near = [beta_near, max(beta_all[i_good,[i_sc,j_sc,k_sc],*], dim=2, /nan)]
						al_compare = [al_compare, al_all[i_good]]
						al_kyoto = [al_kyoto, kal_all[i_good]]
						imf = [imf, imf_all[i_good,*]]
						probes_3p = [probes_3p, rebin_str([[probes[i_sc]], [probes[j_sc]], [probes[k_sc]]], n_good, 3)]

						;;;; compute the distance from the neutral sheet
						l_harris_this = fltarr(n_good)+!values.f_nan
						z_sigma_this = fltarr(n_good)+!values.f_nan
						l_sigma_this = fltarr(n_good)+!values.f_nan
						atanhh = atanh(h_these)
						for i_t = 0, n_good-1 do begin
							;; equation: z = L*atanh(h)+z0
							xls = atanhh[i_t, *] ;; x for least square fit
							yls = pos_these[i_t, *, 2]
							ls_out = poly_fit(xls, yls, 1, sigma = sigma) ;; output is: [z0, L]
							l_harris_this[i_t] = ls_out[1]
							z_sigma_this[i_t] = sigma[0]
							l_sigma_this[i_t] = sigma[1]
						endfor
						l_harris = [l_harris, l_harris_this]
						z_sigma = [z_sigma, z_sigma_this]
						l_sigma = [l_sigma, l_sigma_this]
						zns_these = rebin(l_harris_this, size(h_these, /dim))*atanhh
						zns = [zns, mean(zns_these, dim=2, /nan)]
						zns_max_abs = max(abs(zns_these), i_max_abs, dim=2, /nan)
						zns_far = [zns_far, zns_max_abs*sign(zns_these[i_max_abs])]
						zns_min_abs = min(abs(zns_these), i_min_abs, dim=2, /nan)
						zns_near = [zns_near, zns_min_abs*sign(zns_these[i_min_abs])]
					endif ;; if of having good points
				endfor ;; for of k_sc, the third sc
			endif else begin
				;;;;;;; compute from 2 probes
				;;; first, if 3p is set, coming here means there are only two probes available. So just continue
				if strcmp(method_suf, '_3p') then continue
				tell_good_pos = pos_inrange[*,i_sc] and pos_inrange[*,j_sc] ;; good location
				tell_good_dz = (abs(d_pos[*,2]) lt dsep_max) and (sqrt(d_pos[*,0]^2+d_pos[*,1]^2) lt dother_ratio2sep*abs(d_pos[*,2])) ;; Z separation, good for jx and jy
				tell_good_dy = (abs(d_pos[*,1]) lt dsep_max) and (sqrt(d_pos[*,0]^2+d_pos[*,2]^2) lt dother_ratio2sep*abs(d_pos[*,1])) ;; Y separation, good for jx
				i_good = where(tell_good_pos and (tell_good_dz or tell_good_dy), n_good)

				if n_good gt 0 then begin
					;;;; add to all result arrays
					i_events = [i_events, replicate(i, n_good)]
					t_multi = [t_multi, t_all[i_good]]
					h_i = h_all[i_good,i_sc]
					h_j = h_all[i_good,j_sc]
					hxy_i = hxy_all[i_good,i_sc]
					hxy_j = hxy_all[i_good,j_sc]

					;;;; compute current 
					d_pos_good = d_pos[i_good, *]
					db_good = reform(b_all[i_good, i_sc, *]-b_all[i_good, j_sc, *], n_good, 3)
					b_ave_this = reform(0.5*(b_all[i_good, i_sc, *]+b_all[i_good, j_sc, *]), n_good, 3)
					
					;;;; compute current based on different conditions
					jy_this = fltarr(n_good)+!values.f_nan
					jx_this = fltarr(n_good)+!values.f_nan
					;;; using Y separation
					i_good_dy = where((abs(d_pos_good[*,1]) lt dsep_max) and (sqrt(d_pos_good[*,0]^2+d_pos_good[*,2]^2) lt dother_ratio2sep*abs(d_pos_good[*,1])), n_good_dy)
					if n_good_dy gt 0 then begin
						jx_this[i_good_dy] = db_good[i_good_dy,2]/(mu0*d_pos_good[i_good_dy,1]*RE*1000)
					endif
					;;; using Z separation
					i_good_dz = where((abs(d_pos_good[*,2]) lt dsep_max) and (sqrt(d_pos_good[*,0]^2+d_pos_good[*,1]^2) lt dother_ratio2sep*abs(d_pos_good[*,2])), n_good_dz)
					if n_good_dz gt 0 then begin
						jy_this[i_good_dz] = db_good[i_good_dz,0]/(mu0*d_pos_good[i_good_dz,2]*RE*1000)
						jx_this[i_good_dz] = -db_good[i_good_dz,1]/(mu0*d_pos_good[i_good_dz,2]*RE*1000)
					endif

					jy = [jy, jy_this]
					jx = [jx, jx_this]

					dxyz = [dxyz, abs(d_pos_good)]
					pos_ave = [pos_ave, 0.5*reform(pos_all[i_good, i_sc, *]+pos_all[i_good, j_sc, *], n_good, 3)]
					b_ave = [b_ave, b_ave_this]
					h_these = [[h_i], [h_j]]
					h_ave = [h_ave, 0.5*(h_i+h_j)]
					h_max_abs = max(abs(h_these), i_max_abs, dim=2, /nan)
					h_far = [h_far, h_max_abs*sign(h_these[i_max_abs])]
					h_min_abs = min(abs(h_these), i_min_abs, dim=2, /nan)
					h_near = [h_near, h_min_abs*sign(h_these[i_min_abs])]
					hxy_these = [[hxy_i], [hxy_j]]
					hxy_ave = [hxy_ave, 0.5*(hxy_i+hxy_j)]
					hxy_max_abs = max(abs(hxy_these), i_max_abs, dim=2, /nan)
					hxy_far = [hxy_far, hxy_max_abs*sign(hxy_these[i_max_abs])]
					hxy_min_abs = min(abs(hxy_these), i_min_abs, dim=2, /nan)
					hxy_near = [hxy_near, hxy_min_abs*sign(hxy_these[i_min_abs])]
					beta_ave = [beta_ave, abs(0.5*(beta_all[i_good,i_sc]*sign(h_i)+beta_all[i_good,j_sc]*sign(h_j)))]
					beta_far = [beta_far, min([beta_all[i_good,i_sc],beta_all[i_good,j_sc]])]
					beta_near = [beta_near, max([beta_all[i_good,i_sc],beta_all[i_good,j_sc]])]
					al_kyoto = [al_kyoto, kal_all[i_good]]
					imf = [imf, imf_all[i_good,*]]

					;; rebin the probes array
					probes_dual_this = strarr(n_good, 2)
					for k = 0, n_good-1 do begin
						probes_dual_this[k,*] = transpose(probes[[i_sc,j_sc]])
					endfor ;; for of k, elements of good points
					probes_dual = [probes_dual, probes_dual_this]

					;;;; half thickness of cross-tail current
					Bl_i = Bl[i_good, i_sc]
					Bl_j = Bl[i_good, j_sc]
					Bl_use = fltarr(n_good)

					;; use the Bl of the higher h
					k_use_i = where(abs(h_i) ge abs(h_j), n_use_i)
					if n_use_i gt 0 then Bl_use[k_use_i] = Bl_i[k_use_i]
					k_use_j = where(abs(h_i) lt abs(h_j), n_use_j)
					if n_use_j gt 0 then Bl_use[k_use_j] = Bl_j[k_use_j]
					l_harris_this = d_pos_good[*,2]/(atanh(b_all[i_good, i_sc, 0]/Bl_use)-atanh(b_all[i_good, j_sc, 0]/Bl_use))
					l_harris = [l_harris, l_harris_this]
					Blobe = [Blobe, Bl_use]
					;; compute distance from neutral sheet
					zns_these = l_harris_this*[[atanh(h_i)], [atanh(h_j)]]
					zns = [zns, l_harris_this*0.5*(atanh(h_i)+atanh(h_j))]
					zns_max_abs = max(abs(zns_these), i_max_abs, dim=2, /nan)
					zns_far = [zns_far, zns_max_abs*sign(zns_these[i_max_abs])]
					zns_min_abs = min(abs(zns_these), i_min_abs, dim=2, /nan)
					zns_near = [zns_near, zns_min_abs*sign(zns_these[i_min_abs])]
				endif ;; if of enough good points for 2 probes
			endelse ;; else of 3 or 2 probes
		endfor ;; for of i_sc
	endfor ;; for of j_sc
endfor ;; for of i, events

;;;; trim the result arrays
if n_elements(i_events) gt 1 then i_events = i_events[1:*]
if n_elements(t_multi) gt 1 then t_multi = t_multi[1:*]
if n_elements(pos_ave[*,0]) gt 1 then pos_ave = pos_ave[1:*, *]
if n_elements(b_ave[*,0]) gt 1 then b_ave = b_ave[1:*, *]
if n_elements(beta_ave) gt 1 then beta_ave = beta_ave[1:*]
if n_elements(beta_far) gt 1 then beta_far = beta_far[1:*]
if n_elements(beta_near) gt 1 then beta_near = beta_near[1:*]
if n_elements(h_ave) gt 1 then h_ave = h_ave[1:*]
if n_elements(h_far) gt 1 then h_far = h_far[1:*]
if n_elements(h_near) gt 1 then h_near = h_near[1:*]
if n_elements(hxy_ave) gt 1 then hxy_ave = hxy_ave[1:*]
if n_elements(hxy_far) gt 1 then hxy_far = hxy_far[1:*]
if n_elements(hxy_near) gt 1 then hxy_near = hxy_near[1:*]
if n_elements(al_compare) gt 1 then al_compare = al_compare[1:*]
if n_elements(al_kyoto) gt 1 then al_kyoto = al_kyoto[1:*]
if n_elements(imf[*,0]) gt 1 then imf = imf[1:*, *]
if n_elements(l_harris) gt 1 then l_harris = l_harris[1:*]
if n_elements(zns) gt 1 then zns = zns[1:*]
if n_elements(zns_far) gt 1 then zns_far = zns_far[1:*]
if n_elements(zns_near) gt 1 then zns_near = zns_near[1:*]
if ~strcmp(method_suf, '_3p') then begin
	;; two probes only
	if n_elements(jy) gt 1 then jy = jy[1:*]
	if n_elements(jx) gt 1 then jx = jx[1:*]
	if n_elements(dxyz[*,0]) gt 1 then dxyz = dxyz[1:*, *]
	if n_elements(Blobe) gt 1 then Blobe = Blobe[1:*]
	if n_elements(probes_dual[*,0]) gt 1 then probes_dual = probes_dual[1:*, *]
endif else begin
	;; three probes only
	if n_elements(dir_3p_type) gt 1 then dir_3p_type = dir_3p_type[1:*] ;; 1: close to X; 2: close to Bdir; 3: close to both
	if n_elements(angle_3p_jvx) gt 1 then angle_3p_jvx = angle_3p_jvx[1:*]
	if n_elements(angle_3p_jvb) gt 1 then angle_3p_jvb = angle_3p_jvb[1:*]
	if n_elements(j3p_regx) gt 1 then j3p_regx = j3p_regx[1:*] ;; current value regarding X direction
	if n_elements(j3p_regb) gt 1 then j3p_regb = j3p_regb[1:*] ;; current value regarding B direction
	if n_elements(j3p[*,0]) gt 1 then j3p = j3p[1:*, *]
	if n_elements(probes_3p[*,0]) gt 1 then probes_3p = probes_3p[1:*, *]
	if n_elements(distance_3p[*,0]) gt 1 then distance_3p = distance_3p[1:*, *]
	if n_elements(triangle_3p[*,0]) gt 1 then triangle_3p = triangle_3p[1:*, *]
	if n_elements(z_sigma) gt 1 then z_sigma = z_sigma[1:*]
	if n_elements(l_sigma) gt 1 then l_sigma = l_sigma[1:*]
endelse

;;;; save the data for statistical study
dataout_simple, save_folder+'/i'+list_suf+method_suf, transpose(i_events)
dataout_simple, save_folder+'/t'+list_suf+method_suf, transpose(t_multi)
dataout_simple, save_folder+'/pos'+list_suf+method_suf, transpose(pos_ave)
dataout_simple, save_folder+'/beta'+list_suf+method_suf, transpose(beta_ave)
dataout_simple, save_folder+'/betaFar'+list_suf+method_suf, transpose(beta_far)
dataout_simple, save_folder+'/betaNear'+list_suf+method_suf, transpose(beta_near)
dataout_simple, save_folder+'/h'+list_suf+method_suf, transpose(h_ave)
dataout_simple, save_folder+'/hFar'+list_suf+method_suf, transpose(h_far)
dataout_simple, save_folder+'/hNear'+list_suf+method_suf, transpose(h_near)
dataout_simple, save_folder+'/hxy'+list_suf+method_suf, transpose(hxy_ave)
dataout_simple, save_folder+'/hxyFar'+list_suf+method_suf, transpose(hxy_far)
dataout_simple, save_folder+'/hxyNear'+list_suf+method_suf, transpose(hxy_near)
dataout_simple, save_folder+'/pseudo_al'+list_suf+method_suf, transpose(al_compare)
dataout_simple, save_folder+'/kyoto_al'+list_suf+method_suf, transpose(al_kyoto)
dataout_simple, save_folder+'/b'+list_suf+method_suf, transpose(b_ave)
dataout_simple, save_folder+'/imf'+list_suf+method_suf, transpose(imf)
dataout_simple, save_folder+'/l_harris'+list_suf+method_suf, transpose(l_harris)
dataout_simple, save_folder+'/zns'+list_suf+method_suf, transpose(zns)
dataout_simple, save_folder+'/znsFar'+list_suf+method_suf, transpose(zns_far)
dataout_simple, save_folder+'/znsNear'+list_suf+method_suf, transpose(zns_near)
if ~strcmp(method_suf, '_3p') then begin
	;; two probes only
	dataout_simple, save_folder+'/jy'+list_suf+method_suf, transpose(jy)
	dataout_simple, save_folder+'/jx'+list_suf+method_suf, transpose(jx)
	dataout_simple, save_folder+'/dxyz'+list_suf+method_suf, transpose(dxyz)
	dataout_simple, save_folder+'/Blobe'+list_suf+method_suf, transpose(Blobe)
	output_txt, transpose(probes_dual), filename = save_folder+'/probes_dual'+list_suf+method_suf+'.txt'
endif else begin
	;; three probes only
	dataout_simple, save_folder+'/dir_3p_type'+list_suf, transpose(dir_3p_type)
	dataout_simple, save_folder+'/angle_3p_jvx'+list_suf, transpose(angle_3p_jvx)
	dataout_simple, save_folder+'/angle_3p_jvb'+list_suf, transpose(angle_3p_jvb)
	dataout_simple, save_folder+'/j3p_regx'+list_suf, transpose(j3p_regx)
	dataout_simple, save_folder+'/j3p_regb'+list_suf, transpose(j3p_regb)
	dataout_simple, save_folder+'/j3p'+list_suf, transpose(j3p)
	dataout_simple, save_folder+'/distance_3p'+list_suf, transpose(distance_3p)
	dataout_simple, save_folder+'/triangle_3p'+list_suf, transpose(triangle_3p)
	dataout_simple, save_folder+'/z_sigma'+list_suf+method_suf, transpose(z_sigma)
	dataout_simple, save_folder+'/l_sigma'+list_suf+method_suf, transpose(l_sigma)
	output_txt, transpose(probes_3p), filename = save_folder+'/probes_3p'+list_suf+'.txt'
endelse

;;;; make a list of cases, only time, probe locations judged again later when plotting 
i_sort = sort(t_multi)
t_sort = t_multi[i_sort]

diff = t_sort[1:*]-t_sort[0:-2]
i_break = where(diff gt events_sep*60, n_break)
if n_break gt 1 then begin
	i_break = [i_break, -1]
	events_list = strarr(2, n_break+1)
	;; first
	events_list[*,0] = time_string([t_sort[0], t_sort[i_break[0]]]) 
	;; others
	for i = 1, n_break do begin
		events_list[*,i] = time_string([t_sort[i_break[i-1]+1], t_sort[i_break[i]]]) 
	endfor
endif else begin
	events_list = time_string([t_sort[0], t_sort[-1]])
endelse

;; output the list
output_txt, events_list, filename = save_folder+'/ctcs_events'+list_suf+method_suf+'.txt'

stop
end
