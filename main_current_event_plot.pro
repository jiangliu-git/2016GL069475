pro main_current_event_plot
;;; From the results of main_ctcs_stat, find good events and make plots
thm_init
computer = 'I:'
@folders

pic_folder = pic_folder+'/dual_currents'

;;;; select plot type
;plot_type = 'test' ;; for study
plot_type = 'publication' ;; for paper

;;; select which set of data to use
;list_suf = '' ;; all
list_suf = '_grwexp' ;; late growth (5min) + early expansion phase
;list_suf = '_grwexp20' ;; late growth (20min) + early expansion phase

;;; choose method
;method_suf = ''
method_suf = '_3p' ;; use three probes to compute j in the direction perp to the plane formed by the three

;;; version suffix
;version_suf = '' ;; original
;version_suf = '_v2' ;; new after vassilis comments
version_suf = '_v3' ;; 30 degrees as requirement and nothing else

;;;; choose the theta angle telling jb/jx good or not
theta_max = 30 ;; in degrees, the maximum convincing fac angle (original)
theta_max2 = 45 ;; in degrees, the maximum convincing fac angle (new after vassilis comments)
color_mod = 213 ;; color for cases between 30 and 45

;;;; choose whether to plot jx
;plot_jx = 1 ;; original
plot_jx = 0 ;; new after vassilis comments


;;; set ranges
xrange = [-12, -7.] ;; can show r-1 for growth phase
yrange = [-12, 12.]
zrange = [-10, 10.]
;;; for R-1 current
dawn_end = -2. ;; in RE in Y
dusk_end = 3. ;; in RE in Y

;;; time range of loading (in minutes)
time_load_pre = 20.
time_load_aft = 40.


;;;; load events
case list_suf of
'_sub': listname = 'sb_substorm_list'
else: listname = 'sb'+list_suf+'_list'
endcase
events = load_trange_list(list_folder+'/'+listname+'.txt')

;;;;; load data
i_events = datain_simple(save_folder+'/i'+list_suf+method_suf+'.dat', dim=1, type='int')
t_all = datain_simple(save_folder+'/t'+list_suf+method_suf+'.dat', dim=1, type='double')
pos = datain_simple(save_folder+'/pos'+list_suf+method_suf+'.dat', dim=3, type='float') ;; in RE
h = datain_simple(save_folder+'/h'+list_suf+method_suf+'.dat', dim=1, type='float')

if ~strcmp(method_suf, '_3p') then begin
	;;;; two probes only
	jx = datain_simple(save_folder+'/jx'+list_suf+method_suf+'.dat', dim=1, type='float')
	probes_dual = load_trange_list(save_folder+'/probes_dual'+list_suf+method_suf+'.txt', tempfile = 'sc2_template.sav')
	;!!!!! choose events to examine !!!!!! cannot set i_events_good because the new records are different from the old.
	trange_events = [['2010 4 3 7 58','2010 4 3 8 50'], $ ;ae
;	['2011 5 30 2 10', '2011 5 30 4 10'], $ ;ae
	['2012 7 19 5 38', '2012 7 19 8'], $ ;ae
;	['2012 7 31 4 27', '2012 7 31 6']$ ;de
	['2009 4 11 8 40', '2009 4 11 10 25']$ ;ae
	]
	abc = ['c', 'a', 'b']
	event_num_strs = ['3', '4', '5']
	probes_custom = ['ae', 'ae', 'ae']
endif else begin
	;;; three probes only
	j3p_regx = datain_simple(save_folder+'/j3p_regx'+list_suf+'.dat', dim=1, type='float') ;; current value regarding X direction
	j3p_regb = datain_simple(save_folder+'/j3p_regb'+list_suf+'.dat', dim=1, type='float') ;; current value regarding X direction
	probes_dual = load_trange_list(save_folder+'/probes_3p'+list_suf+'.txt', dim = 3, tempfile = 'sc3_template.sav') ;; though this is 3p, still call it probes_dual for convenience
	;!!!!! choose events to examine !!!!
	i_events_good = [616, 640, 1613]
	;;;; original
	;abc = ['b', 'a', 'c']
	;event_num_strs = ['2', '1', '6']
	;;; new reduced
	n_panels = 7
	event_num_strs = ['2', '1', '3']
	trange_events = [['9 1 31 6 25','9 1 31 7 26'], ['9 2 22 7 55','9 2 22 8 57'], ['12 7 3 2 38','12 7 3 3 47']]
	y_abc = rebin([0.92, 0.8, 0.693, 0.5, 0.467, 0.32, 0.21], n_panels, n_elements(trange_events[0,*]))
endelse


;;; find good events to plot.
;; R-1 current towards earth
tell_inrange = (pos[0,*] gt min(xrange)) and (pos[0,*] lt max(xrange)) and (pos[1,*] gt min(yrange)) and (pos[1,*] lt max(yrange)) and (pos[2,*] gt min(zrange)) and (pos[2,*] lt max(zrange))
tell_dawn = pos[1,*] lt dawn_end
tell_dusk = pos[1,*] gt dusk_end
if ~strcmp(method_suf, '_3p') then begin
	tell_r1_dawn = (jx gt 0) and tell_dawn ;; towards earth
	tell_r1_dusk = (jx lt 0) and tell_dusk ;; away from earth
endif else begin
	tell_r1_dawn = (j3p_regx gt 0) and tell_dawn ;; towards earth
	tell_r1_dusk = (j3p_regx lt 0) and tell_dusk ;; away from earth
endelse
store_data, '_r1dawn', data = tell_r1_dawn
store_data, '_r1dusk', data = tell_r1_dusk
dirs = ['_r1dawn', '_r1dusk']

i_epub = 0 ;; for publication only
for i = 0, n_elements(dirs)-1 do begin
	get_data, dirs[i], data = tell_loc
	k_loc = where(tell_loc and tell_inrange, n_loc)
	if n_loc lt 1 then continue
	
	i_events_loc_all = i_events[k_loc]
	i_events_loc = i_events_loc_all[uniq(i_events_loc_all, sort(i_events_loc_all))]

	for i_event = 0, n_elements(i_events_loc)-1 do begin
		trange = time_double(events[*,i_events_loc[i_event]])
		trange_show = trange+[-time_load_pre, time_load_aft]*60. ;; range to load and plot the data

		;;; using tranges instead of event number, works for both 3p and 2p
		if keyword_set(trange_events) then begin
			;;; see whether there is a match
			ifmatch = 0
			for i_trange_cust = 0, n_elements(trange_events[0,*])-1 do begin
				if (~strcmp(trange_events[0, i_trange_cust], '')) and (~strcmp(trange_events[1, i_trange_cust], '')) then begin
					trange_this = time_double(trange_events[*,i_trange_cust])
					if ((trange_this[0] ge trange_show[0]) and (trange_this[0] le trange_show[1])) or ((trange_this[1] ge trange_show[0]) and (trange_this[1] le trange_show[1])) or ((trange_show[0] ge trange_this[0]) and (trange_show[0] le trange_this[1])) or ((trange_show[1] ge trange_this[0]) and (trange_show[1] le trange_this[1])) then begin
						trange_show = trange_this
						i_epub = i_trange_cust
						ifmatch = 1
						break
					endif
				endif
			endfor
			if ~ifmatch then continue
		endif

		if keyword_set(i_events_good) then begin
			;;; plot only good events specified manually
			if total(i_events_loc[i_event] eq i_events_good) eq 0 then continue
			if keyword_set(trange_events) then begin
				if ~strcmp(trange_events[0, i_epub], '') then trange_show[0] = time_double(trange_events[0, i_epub])
				if ~strcmp(trange_events[1, i_epub], '') then trange_show[1] = time_double(trange_events[1, i_epub])
			endif
		endif

		trange_load = trange_show+[-100., 100.]

		del_data, 'th*'
		;; load al and initiate the tplot variables
		load_bin_data, datatype = 'kyoto_al', trange = trange_load, datafolder = kyoto_al_folder, /tclip
		options, 'kyoto_al_tclip', ytitle = 'AL [nT]', ysubtitle = ''
		tnames = 'kyoto_al_tclip'
		;; find the probes and range of good j
		k_this = where(i_events eq i_events_loc[i_event], n_this)
		if n_this lt 1 then message, 'Impossible!'
		probes_dual_this = probes_dual[*,k_this]
		if strcmp(method_suf, '') and keyword_set(probes_custom) then begin
			probes_plot = [strmid(probes_custom[i_epub],0,1), strmid(probes_custom[i_epub],1,1)]
		endif else begin
			probes_plot = probes_dual_this[uniq(probes_dual_this, sort(probes_dual_this))]
		endelse
		;; load position and magnetic field and add to tplot variables
		load_bin_data, probes = probes_plot, datatype = 'pos', trange = trange_load, datafolder = pos_folder, /tclip
		load_bin_data, probes = probes_plot, datatype = 'fgs', trange = trange_load, datafolder = fgs_folder
		load_bin_data, probes = probes_plot, datatype = 'Pall', trange = trange_load, datafolder = Pall_folder, /tclip
		;; load newest magnetic field
		;thm_load_fit, probes = probes_plot, datatype = 'fgs', trange = trange_load, coord='gsm', level = 1, suffix = '_gsm'
		for i_p = 0, n_elements(probes_plot)-1 do begin
			;;; adjust fgs data for perticular events
			case probes_plot[i_p] of
				'a': bad_trange = ['9 1 31 7 21 25', '9 1 31 7 21 45']
				'd': bad_trange = ['9 1 31 7 15 45', '9 1 31 7 16 20']
				else: bad_trange = ''
			endcase
			if ~strcmp(bad_trange[0], '') then begin
				get_data, 'th'+probes_plot[i_p]+'_fgs_gsm', t_adjust, fgs_adjust
				i_bad = where((t_adjust gt time_double(bad_trange[0])) and (t_adjust lt time_double(bad_trange[1])), n_bad)
				if n_bad gt 0 then begin
					for i_row = 0, 2 do begin
						t_new = [t_adjust[0:i_bad[0]-1], t_adjust[i_bad[-1]+1:*]]
						fgs_new = [fgs_adjust[0:i_bad[0]-1, *], fgs_adjust[i_bad[-1]+1:*, *]]
					endfor
					store_data, 'th'+probes_plot[i_p]+'_fgs_gsm', data = {x:t_new, y:fgs_new}
				endif
			endif
			time_clip, 'th'+probes_plot[i_p]+'_fgs_gsm', trange_load[0], trange_load[1]
			options, 'th'+probes_plot[i_p]+'_fgs_gsm*', ytitle = 'B!d'+thm_probe_color(probes_plot[i_p], /num)+'!n [nT]', ysubtitle = '', colors = [2,4,6], labels = ['B!dx', 'B!dy', 'B!dz'], labflag = 1
		endfor
		tnames = [tnames, 'th'+probes_plot+'_fgs_gsm_tclip']
		;; compute scale height
		;; find out the dual computations and add j to tplot variables, and the range of good j and mark them
		if strcmp(method_suf, '_3p') then begin
			probes_comb = probes_dual_this[0,*]+probes_dual_this[1,*]+probes_dual_this[2,*]
		endif else begin
			probes_comb = probes_dual_this[0,*]+probes_dual_this[1,*]
		endelse
		if strcmp(plot_type, 'test') then begin
			if strcmp(method_suf, '_3p') then begin
				tnames_bar = ['','','','','','']
			endif else begin
				tnames_bar = ['','','','']
			endelse
		endif else tnames_bar = ''
		if keyword_set(probes_custom) then begin
			probes_uniq = probes_custom[i_epub]
		endif else begin
			probes_uniq = probes_comb[uniq(probes_comb, sort(probes_comb))]
		endelse
		t_this = t_all[k_this]
		varlabels = ''
		for i_probes = 0, n_elements(probes_uniq)-1 do begin
			;;;;;; compute j
			probes_comb_this = probes_uniq[i_probes]
			sc1 = strmid(probes_comb_this, 0, 1)
			sc2 = strmid(probes_comb_this, 1, 1)
			;;;; all interpol to the location of sc1
			;interpol_tname = 'th'+sc1+'_fgs_gsm_tclip' ;; 3s resolution
			interpol_tname =  'th'+sc1+'_state_pos_tclip' ;; 1min resolution
			tnames_interpol = ['th'+sc1+'_state_pos_tclip', 'th'+sc2+'_state_pos_tclip', 'th'+sc1+'_fgs_gsm_tclip', 'th'+sc2+'_fgs_gsm_tclip', 'th'+sc1+'_Pall_tclip', 'th'+sc2+'_Pall_tclip']
			for i_name = 0, n_elements(tnames_interpol)-1 do begin
				tinterpol_mxn, tnames_interpol[i_name], interpol_tname
				;;;;; Mark all interpolated data's NaN points (where interpolated fake points exist)
				get_data, tnames_interpol[i_name]+'_interp', t, data_this
				get_data, tnames_interpol[i_name], t_orig, nouse
				i_out = where((t lt t_orig[0]) or (t gt t_orig[-1]), n_out)
				if n_out gt 0 then data_this[i_out, *] = !values.f_nan
				store_data, tnames_interpol[i_name]+'_interp', data = {x:t, y:data_this}
			endfor
			get_data, 'th'+sc1+'_state_pos_tclip_interp', t, pos1
			get_data, 'th'+sc2+'_state_pos_tclip_interp', t, pos2
			get_data, 'th'+sc1+'_fgs_gsm_tclip_interp', t, b1
			get_data, 'th'+sc2+'_fgs_gsm_tclip_interp', t, b2
			get_data, 'th'+sc1+'_Pall_tclip_interp', t, Pall1
			get_data, 'th'+sc2+'_Pall_tclip_interp', t, Pall2
			;;;;; compute scale height: proxy distance from neutral sheet
			if keyword_set(Pall1) then h1 = b1[*,0]/sqrt((b1[*,0]^2+b1[*,1]^2)+Pall1[*,1]/nTesla2_to_nPa) else h1 = !values.f_nan
			if keyword_set(Pall2) then h2 = b2[*,0]/sqrt((b2[*,0]^2+b2[*,1]^2)+Pall2[*,1]/nTesla2_to_nPa) else h2 = !values.f_nan
			;;; compute currents and locations
			if strcmp(method_suf, '_3p') then begin
				;;; load the 3rd spacecraft data
				sc3 = strmid(probes_comb_this, 2, 1)
				tinterpol_mxn,'th'+sc3+'_state_pos_tclip', interpol_tname
				tinterpol_mxn,'th'+sc3+'_fgs_gsm_tclip', interpol_tname
				tinterpol_mxn,'th'+sc3+'_Pall_tclip', interpol_tname
				get_data, 'th'+sc3+'_state_pos_tclip_interp', t, pos3
				get_data, 'th'+sc3+'_fgs_gsm_tclip_interp', t, b3
				get_data, 'th'+sc3+'_Pall_tclip_interp', t, Pall3
				;;;;; compute scale height: proxy distance from neutral sheet
				if keyword_set(Pall3) then h3 = b3[*,0]/sqrt((b3[*,0]^2+b3[*,1]^2)+Pall3[*,1]/nTesla2_to_nPa) else h3 = !values.f_nan
				;;;;; compute current and locations
				;;; compute the current vector no matter what, this is a vector, in the direction of the plane formed by the three probes (make it 3xn_t_arr)
				j_this = curlometer3p(transpose(b1), transpose(b2), transpose(b3), $
							transpose(pos1), transpose(pos2), transpose(pos3), $
							distance_12 = distance_12, distance_23 = distance_23, distance_31 = distance_31, position_mean = position_mean, $
							angle_12 = angle_12, angle_23 = angle_23, angle_31 = angle_31)

				b_ave = (b1+b2+b3)/3 ;; the average direction of the field lines of the three spacecraft
				angle_j_x = angle_vectors(j_this, rebin([1.,0,0], 3, n_elements(t)), /deg, /ignore_dir)
				angle_j_b = angle_vectors(j_this, transpose(b_ave), /deg, /ignore_dir)
				j_str = sqrt(total(j_this^2, 1))
				sign_j_x = transpose(sign(dotp_long(j_this, rebin([1.,0,0], 3, n_elements(t)))))
				sign_j_b = transpose(sign(dotp_long(j_this, transpose(b_ave))))*sign(b_ave[*,0])
				j3p_regx = j_str*sign_j_x
				j3p_regb = j_str*sign_j_b
				;;; average locations
				pos_ave_re = (pos1+pos2+pos3)/(3*RE)
				h_ave = (h1+h2+h3)/3
			endif else begin
				;;;;; compute current and locations
				dpos = pos2-pos1
				pos_ave_re = 0.5*(pos2+pos1)/RE
				db = b2-b1
				jx = -db[*,1]/(mu0*dpos[*,2]*1000)
				jy = db[*,0]/(mu0*dpos[*,2]*1000)
				;;; average scale height
				h_ave = 0.5*(h1+h2)
			endelse
			;;;;; save and edit data
			if strcmp(method_suf, '_3p') then begin
				distances = transpose([distance_12, distance_23, distance_31])/RE
				triangles = transpose([angle_12, angle_23, angle_31])
				store_data, 'th'+probes_comb_this+'_seperations', data = {x:t, y:distances}
				store_data, 'th'+probes_comb_this+'_sepmax', data = {x:t, y:max(distances, dim=2)}
				store_data, 'th'+probes_comb_this+'_triangles', data = {x:t, y:triangles}
				store_data, 'th'+probes_comb_this+'_anglemin', data = {x:t, y:min(triangles, dim=2)}
				if plot_jx then begin
					store_data, 'th'+probes_comb_this+'_j_angles', data = {x:t, y:transpose([angle_j_x, angle_j_b])}
					options, 'th'+probes_comb_this+'_j_angles', ytitle = theta_letter+'!dj!n [degs]', colors = [2, 1], labels = theta_letter+'!dj'+['X', 'B']+'!n', labflag = 1
				endif else begin
					store_data, 'th'+probes_comb_this+'_j_angles', data = {x:t, y:transpose(angle_j_b)}
					options, 'th'+probes_comb_this+'_j_angles', ytitle = theta_letter+'!djB!n [degs]'
				endelse
				if strcmp(plot_type, 'publication') then begin
					;; manipulate reg B
					j3p_regbr = -j3p_regb*sign(pos_ave_re[*,1])
					j3p_regbr_all = j3p_regbr
					j3p_regbr_mod = j3p_regbr
					j3p_regbr_bad = j3p_regbr
					i_bad_regb = where(angle_j_b gt theta_max, n_bad_regb)
					i_mod_regb = where(angle_j_b gt theta_max2, n_mod_regb)
					i_good_regb = where(angle_j_b le theta_max, n_good_regb)
					if n_bad_regb gt 0 then j3p_regbr[i_bad_regb] = !values.f_nan
					if n_mod_regb gt 0 then j3p_regbr_mod[i_mod_regb] = !values.f_nan
					if n_good_regb gt 0 then j3p_regbr_bad[i_good_regb] = !values.f_nan
					;; manipulate reg X
					j3p_regxr = -j3p_regx*sign(pos_ave_re[*,1])
					j3p_regxr_all = j3p_regxr  
					j3p_regxr_bad = j3p_regxr  
					i_bad_regx = where(angle_j_x gt theta_max, n_bad_regx)
					i_good_regx = where(angle_j_x le theta_max, n_good_regx)
					if n_bad_regx gt 0 then j3p_regxr[i_bad_regx] = !values.f_nan
					if n_good_regx gt 0 then j3p_regxr_bad[i_good_regx] = !values.f_nan
					;; store data
					if plot_jx then begin
						store_data, 'th'+probes_comb_this+'_j_good', data = {x:t, y:[[j3p_regxr], [j3p_regbr]]}
						store_data, 'th'+probes_comb_this+'_j_bad', data = {x:t, y:[[j3p_regxr_bad], [j3p_regbr_bad]]}
					endif else begin
						store_data, 'th'+probes_comb_this+'_j_good', data = {x:t, y:j3p_regbr}
						store_data, 'th'+probes_comb_this+'_j_mod', data = {x:t, y:j3p_regbr_mod}
						;store_data, 'th'+probes_comb_this+'_j_bad', data = {x:t, y:j3p_regbr_bad}
						store_data, 'th'+probes_comb_this+'_j_bad', data = {x:t, y:j3p_regbr_all} ;; use all to get well -connected plots
					endelse
					if strcmp(version_suf, '_v2') then begin
						store_data, 'th'+probes_comb_this+'_j', data = 'th'+probes_comb_this+'_j_'+['bad', 'mod', 'good']
					endif else begin
						store_data, 'th'+probes_comb_this+'_j', data = 'th'+probes_comb_this+'_j_'+['bad', 'good']
					endelse
					options, 'th'+probes_comb_this+'_j_good', thick = 1.1
					options, 'th'+probes_comb_this+'_j_mod', color = color_mod
					options, 'th'+probes_comb_this+'_j_bad', linestyle = 1
					if plot_jx then begin
						options, 'th'+probes_comb_this+'_j_good', colors = [3, 1], labels = ['j!s!dx!i*!r!uR!n', 'j!s!d//!r!uR!n'], labflag = 1
						options, 'th'+probes_comb_this+'_j_bad', colors = [3, 1]
						options, 'th'+probes_comb_this+'_j', ytitle = 'j!uR!n [nA/m!u2!n]', ysubtitle = ''
					endif else begin
						options, 'th'+probes_comb_this+'_j', ytitle = 'j!s!d//!r!uR!n [nA/m!u2!n]', ysubtitle = ''
					endelse
				endif else begin
					store_data, 'th'+probes_comb_this+'_j', data = {x:t, y:[[j3p_regx], [j3p_regb]]}
					options, 'th'+probes_comb_this+'_j', ysubtitle = '[nA/m!u2!n]', colors = [2, 1], labels = ['regX', 'regB'], labflag = 1
				endelse
				options, 'th'+probes_comb_this+'_seperations', ysubtitle = '[RE]', colors = [2, 4, 6], labels = [sc1+sc2, sc2+sc3, sc3+sc1], labflag = 1
				options, 'th'+probes_comb_this+'_triangles', ysubtitle = '[deg]', colors = [2, 4, 6], labels = [sc1+sc2, sc2+sc3, sc3+sc1], labflag = 1
				ylim, 'th'+probes_comb_this+'_j_angles', 0, 90.
				options, 'th'+probes_comb_this+'_sepmax', ytitle = 'D!dmax!n [R!dE!n]'
				options, 'th'+probes_comb_this+'_anglemin', ytitle = alpha_letter+'!dmin!n [degs]'
			endif else begin
				store_data, 'th'+probes_comb_this+'_dpos_re', data = {x:t, y:[[sqrt(total(dpos[*,0:1]^2, 2))], [abs(dpos[*,2])]]/RE}
				if strcmp(plot_type, 'publication') then begin
					jxr = -jx*sign(pos_ave_re[*,1])
					store_data, 'th'+probes_comb_this+'_j', data = {x:t, y:jxr}
					options, 'th'+probes_comb_this+'_dpos_re', ytitle = 'd [R!dE!n]', ysubtitle = '', colors = [78, 6], labels = ['dXY', 'dZ'], labflag = 1
					options, 'th'+probes_comb_this+'_j', ytitle = 'j!s!dx!r!uR!n [nA/m!u2!n]', colors = 2
				endif else begin
					store_data, 'th'+probes_comb_this+'_j', data = {x:t, y:[[jx], [jy]]}
					options, 'th'+probes_comb_this+'_dpos_re', ysubtitle = '[R!dE!n]', colors = [6, 1], labels = ['dZ', 'dXY']
					options, 'th'+probes_comb_this+'_j', ysubtitle = '[nA/m!u2!n]', colors = [2, 4], labels = ['j!dx', 'j!dy'], labflag = 1
				endelse
			endelse
			store_data, 'th'+probes_comb_this+'_x_ave_re', data = {x:t, y:pos_ave_re[*,0]}
			store_data, 'th'+probes_comb_this+'_y_ave_re', data = {x:t, y:pos_ave_re[*,1]}
			store_data, 'th'+probes_comb_this+'_h_ave', data = {x:t, y:h_ave}
			if strcmp(plot_type, 'publication') then begin
				options, 'th'+probes_comb_this+'_x_ave_re', ytitle = jiao_l+'X'+jiao_r+'!n [R!dE!n]', ysubtitle = ''
				options, 'th'+probes_comb_this+'_y_ave_re', ytitle = jiao_l+'Y'+jiao_r+'!n [R!dE!n]', ysubtitle = ''
				options, 'th'+probes_comb_this+'_h_ave', ytitle = jiao_l+'h'+jiao_r, ysubtitle = ''
			endif else begin
				options, 'th'+probes_comb_this+'_x_ave_re', ytitle = 'X!d'+probes_comb_this+'!n [R!dE!n]', ysubtitle = ''
				options, 'th'+probes_comb_this+'_y_ave_re', ytitle = 'Y'+probes_comb_this+'!n [R!dE!n]', ysubtitle = ''
				options, 'th'+probes_comb_this+'_h_ave', ytitle = 'h!d'+probes_comb_this, ysubtitle = ''
			endelse
			;; tell whether location within range, store if yes
			if strmatch(dirs[i], '*dawn*') then begin
				nouse = where(pos_ave_re[*,1] lt dawn_end, n_within)
				jx_look = 'positive'
			endif
			if strmatch(dirs[i], '*dusk*') then begin
				nouse = where(pos_ave_re[*,1] gt dusk_end, n_within)
				jx_look = 'negative'
			endif
			if n_within gt 0 then begin
				;;; find the labels to write in the figure's bottom
				varlabels = [varlabels, 'th'+probes_comb_this+'_y_ave_re', 'th'+probes_comb_this+'_x_ave_re']
				if strcmp(plot_type, 'publication') then begin
					if strcmp(method_suf, '_3p') then begin
						varlabels = [varlabels, 'th'+probes_comb_this+'_sepmax', 'th'+probes_comb_this+'_anglemin']
					endif
				endif else begin
					if finite(mean(h_ave, /nan)) then varlabels = [varlabels, 'th'+probes_comb_this+'_h_ave']
				endelse
				;; find out the time to mark
				k_comb_this = where(strcmp(probes_comb, probes_uniq[i_probes]), n_comb_this) ;; where in the data record the probes match, meaning good
				if n_comb_this lt 1 then message, 'Impossible!'
				t_comb_this = t_this[k_comb_this] ;; time points when these two satellites are seperated within the range defined in main_ctcs_stat, when j is computed.
				;; tplot names to plot, tplot names to add bars, and times to add the bars of good time
				tnames = [tnames, 'th'+probes_comb_this+'_j']
				if strcmp(plot_type, 'publication') then begin
					tnames = [tnames, 'th'+probes_comb_this+'_h_ave']
					if strcmp(method_suf, '_3p') then begin
						tnames = [tnames, 'th'+probes_comb_this+'_j_angles']
					endif else begin
						tnames = [tnames, 'th'+probes_comb_this+'_dpos_re']
					endelse
				endif else begin
					if strcmp(method_suf, '_3p') then begin
						tnames = [tnames, 'th'+probes_comb_this+'_j_angles', 'th'+probes_comb_this+'_seperations', 'th'+probes_comb_this+'_triangles']
						tnames_bar_this = [tnames_bar_this, 'th'+probes_comb_this+'_j_angles', 'th'+probes_comb_this+'_seperations', 'th'+probes_comb_this+'_triangles']
					endif else begin
						tnames = [tnames, 'th'+probes_comb_this+'_dpos_re']
						tnames_bar_this = [tnames_bar_this, 'th'+probes_comb_this+'_dpos_re']
					endelse
					tnames_bar_this = [tnames_bar_this, time_string([t_comb_this[0], t_comb_this[-1]])]
					tnames_bar = [[tnames_bar], [tnames_bar_this]]
				endelse
			endif
		endfor ;; for of i_probes, the dual observations of each growth/other phase events

		;;; make the plot
		save_path = pic_folder+'/'+strmid(dirs[i], 1)+'_event'+strcompress(string(fix(i_events_loc[i_event])), /remove)+method_suf +version_suf
		if strcmp(plot_type, 'publication') then begin
			popen, save_path
			print_options, xsize=5.5, ysize=10
		endif

		if n_elements(varlabels) gt 0 then varlabels = varlabels[1:*]
		if strcmp(plot_type, 'test') then begin
			title = strmid(dirs[i], 1)+', look for '+jx_look+' j!dx!n. Event'+strcompress(string(fix(i_events_loc[i_event])))
		endif else begin
			title = 'Event '+event_num_strs[i_epub]
		endelse

		tplot, tnames, var_label = varlabels, title = title, trange = trange_show
		timebar_mass, 0, varname = tnames, /databar, line = 2
		if strcmp(plot_type, 'test') then timebar, trange, varname = 'kyoto_al_tclip', line = 1

		;;;; bars or arrows of time
		if strcmp(plot_type, 'test') and (n_elements(tnames_bar[0,*]) gt 1) then begin
			for i_bar = 1, n_elements(tnames_bar[0,*])-1 do begin
				timebar_mass, [tnames_bar[-2,i_bar], tnames_bar[-1,i_bar]], varname = tnames_bar[0:n_elements(tnames_bar[*,0])-3,i_bar], line = 1
			endfor
		endif else begin
			;;; time bars
			timebar, '9 2 22 8 20', line = 1 ;; event 1
		endelse

		;;;; bars of values
		if strcmp(method_suf, '_3p') then begin
			timebar_mass, [theta_max, 180.-theta_max], varname = 'th???_j_angles', /databar, line = 1
			if keyword_set(theta_max2) and strcmp(version_suf, '_v2') then timebar_mass, [theta_max2, 180.-theta_max2], varname = 'th???_j_angles', /databar, line = 1, color = color_mod
			if strcmp(plot_type, 'test') then begin
				timebar_mass, 1.5, varname = 'th???_seperations', /databar, line = 1
				timebar_mass, 15., varname = 'th???_triangles', /databar, line = 1
			endif
		endif

		if strcmp(plot_type, 'publication') then begin
			;;;;; draw arrows
			case i_events_loc[i_event] of
				640: begin ;; event 1
					add_pointer, 0.4, 0.42, '1', dir = 'southeast'
					add_pointer, 0.5, 0.455, '2', dir = 'southeast'
					add_pointer, 0.725, 0.42, '3', dir = 'southwest'
					add_pointer, 0.733, 0.398, '4', dir = 'northeast'
					end
				616: begin ;; event 2
					add_pointer, 0.534, 0.418, '5', dir = 'southeast'
					add_pointer, 0.58, 0.47, '6', dir = 'east'
					add_pointer, 0.629, 0.444, '7', dir = 'west'
					end
				928: begin ;; event 3
					add_pointer, 0.33, 0.488, '8', dir = 'southeast'
					add_pointer, 0.426, 0.348, '9', dir = 'northeast'
					end
				1641: begin ;; event 4
					add_pointer, 0.325, 0.48, '10', dir = 'southeast'
					add_pointer, 0.45, 0.478, '11', dir = 'southwest'
					end

				else: print, 'No arrow drawn for this event'
			endcase

			;;;;; write abcs
			xyouts, replicate(0.22, n_panels), y_abc[*,i_epub], '('+getletters(indgen(n_panels)+n_panels*(float(event_num_strs[i_epub])-1))+')', /normal

			pclose
		endif else begin
			makepng, save_path
		endelse
		
		if keyword_set(i_events_good) then i_epub = i_epub+1
	endfor ;; for of i_event, growth/other phase events
endfor

stop
end
