pro main_current_stat_plot
;;; plot the results from main_ctcs_stat or main_ctcs_stat_lists
thm_init
@folders

;;;;;;;;; choose plot type, test or publication ;;;;
;plot_type = 'test' ;; plots for study
plot_type = 'publication' ;; plots for publication
;plot_type = 'proposal' ;; plots for proposal

;;; select which set of data to use
list_suf = '' ;; all
;list_suf = '_grw_csv' ;; growth phase
;list_suf = '_grwexp' ;; late growth + early expansion phase
;list_suf = '_grwexp20' ;; late growth + early expansion phase
;list_suf = '_nongrwexp' ;; everything but late growth + early expansion phase
;list_suf = '_sub' ;; substorm
;list_suf = '_exp_csv' ;; conservative (early) expansion phase
;list_suf = '_rcv_lbl' ;; liberal recovery phase

;;; select which method of computing data
;method_suf = '' ;; use seperated probes in Z
method_suf = '_3p' ;; use three probes to compute j in the direction perp to the plane formed by the three

;;; select the direction suffix (only work for 3p method)
if strcmp(method_suf, '') then dir_suf = '' else begin
;	dir_suf = '' ;; regarding x
	dir_suf = '_regb' ;; regarding B
;	dir_suf = '_regbstrict' ;; regarding B, but strict: the final j//x and original jx should not have different signs.
endelse

;;; select the jvb angle suffix (default is 30 degrees)
jvbangle_suf = '' ;; 30 degrees
;jvbangle_suf = '_10' ;; 10 degrees
;jvbangle_suf = '_15' ;; 15 degrees
;jvbangle_suf = '_45' ;; 45 degrees
;jvbangle_suf = '_90' ;; 90 degrees, all

;;; set ranges
xrange = [-12, -8.] ;; can show r-1 for growth phase
yrange = [-12, 12.]
zrange = [-10, 10.]
jy_range = [0., 3.] ;;for plot
;jx_range = [-3, 3.] ;;for plot
jx_range = [-0.3, 0.3] ;;for plot
al_low = -50.
if strcmp(plot_type, 'proposal') then begin
	y_dawn = -1.5
	y_dusk = 1.5
endif else begin
	y_dawn = -2
	y_dusk = 2
endelse

;;; more strict requirements for 3p computed currents.
case jvbangle_suf of 
	'_10': angle_3p_jvb_req = 10. ;; more strict value
	'_15': angle_3p_jvb_req = 15. ;; more strict value
	'_45': angle_3p_jvb_req = 45. ;; VA suggested value
	'_90': angle_3p_jvb_req = 90. ;; disable this criterion
	else: angle_3p_jvb_req = 30. ;; in degrees, must be smaller than this, main_current_stat requires 30. Requiring 15 degs will reduce events by 96%
endcase
if strcmp(plot_type, 'proposal') then begin
	triangle_3p_req = 15. 
endif else begin
	triangle_3p_req = 30. ;; in degrees, must be larger than this, main_current_stat requires 15. Requiring 30 degs will reduce events by 65%, but gives better results. Use 30 degs for publication.
endelse
distance_3p_req = 1.5 ;; in RE, must be smaller than this, main_current_stat requires 2. Requiring 1.5 will reduce events by ~9% (use for publication). Requiring 1 will reduce events further by 86%.

;; plot specifics
al_high = -100.
c_pts = 5 ;; minumum points to be plotted
c_events = 5 ;; minumum number of events in each bin to be plotted
binsizey = 6.
binsizez = 1.

;;; load data
t_all = datain_simple(save_folder+'/t'+list_suf+method_suf+'.dat', dim=1, type='double')
i_all = datain_simple(save_folder+'/i'+list_suf+method_suf+'.dat', dim=1, type='int') ;; note: this event is just the splitted time ranges, not real event.
pos = datain_simple(save_folder+'/pos'+list_suf+method_suf+'.dat', dim=3, type='float')
beta_ave = datain_simple(save_folder+'/beta'+list_suf+method_suf+'.dat', dim=1, type='float')
beta_far = datain_simple(save_folder+'/betaFar'+list_suf+method_suf+'.dat', dim=1, type='float')
beta_near = datain_simple(save_folder+'/betaNear'+list_suf+method_suf+'.dat', dim=1, type='float')
b_ave = datain_simple(save_folder+'/b'+list_suf+method_suf+'.dat', dim=3, type='float')
h = datain_simple(save_folder+'/h'+list_suf+method_suf+'.dat', dim=1, type='float')
h_far = datain_simple(save_folder+'/hFar'+list_suf+method_suf+'.dat', dim=1, type='float')
h_near = datain_simple(save_folder+'/hNear'+list_suf+method_suf+'.dat', dim=1, type='float')
hxy = datain_simple(save_folder+'/hxy'+list_suf+method_suf+'.dat', dim=1, type='float')
hxy_far = datain_simple(save_folder+'/hxyFar'+list_suf+method_suf+'.dat', dim=1, type='float')
hxy_near = datain_simple(save_folder+'/hxyNear'+list_suf+method_suf+'.dat', dim=1, type='float')
;al = datain_simple(save_folder+'/pseudo_al'+list_suf+method_suf+'.dat', dim=1, type='double') ;; original
al = datain_simple(save_folder+'/kyoto_al'+list_suf+method_suf+'.dat', dim=1, type='double')
imf = datain_simple(save_folder+'/imf'+list_suf+method_suf+'.dat', dim=3, type='float')
l_harris = datain_simple(save_folder+'/l_harris'+list_suf+method_suf+'.dat', dim=1, type='float') ;; half width of harris sheet
zns = datain_simple(save_folder+'/zns'+list_suf+method_suf+'.dat', dim=1, type='float') ;; distance from the neutral sheet in RE
zns_far = datain_simple(save_folder+'/znsFar'+list_suf+method_suf+'.dat', dim=1, type='float') ;; distance from the neutral sheet in RE
zns_near = datain_simple(save_folder+'/znsNear'+list_suf+method_suf+'.dat', dim=1, type='float') ;; distance from the neutral sheet in RE
rho = sqrt(pos[0,*]^2+pos[1,*]^2) ;; rho = sqrt(x2+y2)
if ~strcmp(method_suf, '_3p') then begin
	;;;; two probes only
	dxyz = datain_simple(save_folder+'/dxyz'+list_suf+method_suf+'.dat', dim=3, type='float') ;; using this can further limit the seperation range
	jx = datain_simple(save_folder+'/jx'+list_suf+method_suf+'.dat', dim=1, type='float')
	jy = datain_simple(save_folder+'/jy'+list_suf+method_suf+'.dat', dim=1, type='float')
	Blobe = datain_simple(save_folder+'/Blobe'+list_suf+method_suf+'.dat', dim=1, type='float') ;; half width of harris sheet
	probes_dual = load_trange_list(save_folder+'/probes_dual'+list_suf+method_suf+'.txt', tempfile = 'sc2_template.sav')
	;; assign the fac current
	jfac = jx
endif else begin
	;;; three probes only
	j3p_regx = datain_simple(save_folder+'/j3p_regx'+list_suf+'.dat', dim=1, type='float') ;; current value regarding X direction
	j3p_regb = datain_simple(save_folder+'/j3p_regb'+list_suf+'.dat', dim=1, type='float') ;; current value regarding X direction
	j3p = datain_simple(save_folder+'/j3p'+list_suf+'.dat', dim=3, type='float') ;; current value regarding X direction
	distance_3p = datain_simple(save_folder+'/distance_3p'+list_suf+'.dat', dim=3, type='float') ;; satellite seperations
	triangle_3p = datain_simple(save_folder+'/triangle_3p'+list_suf+'.dat', dim=3, type='float') ;; angles of the triangle formed by the 3 probes.
	angle_3p_jvx = datain_simple(save_folder+'/angle_3p_jvx'+list_suf+'.dat', dim=1, type='float') ;; angles of j and X
	angle_3p_jvb = datain_simple(save_folder+'/angle_3p_jvb'+list_suf+'.dat', dim=1, type='float') ;; angles of j and B
	z_sigma = datain_simple(save_folder+'/z_sigma'+list_suf+method_suf+'.dat', dim=1, type='float') ;; half width of harris sheet
	l_sigma = datain_simple(save_folder+'/l_sigma'+list_suf+method_suf+'.dat', dim=1, type='float') ;; half width of harris sheet
	probes_dual = load_trange_list(save_folder+'/probes_3p'+list_suf+'.txt', dim = 3, tempfile = 'sc3_template.sav') ;; though this is 3p, still call it probes_dual for convenience
	;; note: if away from a certain vector it will be NaN
	case dir_suf of
	'': begin
		jfac = j3p_regx
		angle_vec = angle_3p_jvx
		end
	'_regb': begin
		jfac = j3p_regb
		angle_vec = angle_3p_jvb
		end
	'_regbstrict': begin
		jfac = j3p_regb
		angle_vec = angle_3p_jvb
		end
	endcase
endelse

;;;;;;;; locations or conditions
tell_all = al lt 1 ;; all points
tell_al_low = al gt al_low
tell_al_high = al lt al_high
tell_dawn = pos[1,*] lt y_dawn
tell_dusk = pos[1,*] gt y_dusk
tell_south = pos[2,*] lt 0
tell_north = pos[2,*] gt 0
tell_imfsouth = imf[2,*] lt 0.
tell_imfnorth = imf[2,*] gt 0.
tell_inx = (pos[0,*] gt min(xrange)) and (pos[0,*] lt max(xrange))
tell_iny = (pos[1,*] gt min(yrange)) and (pos[1,*] lt max(yrange))
tell_inz = (pos[2,*] gt min(zrange)) and (pos[2,*] lt max(zrange))
tell_ncps = abs(b_ave[0,*]) gt abs(b_ave[2,*]) ;; not central plasma sheet: |Bx|>|Bz|
if strcmp(method_suf, '_3p') then begin
	tell_good_sep3 = (distance_3p[0,*] lt distance_3p_req) and (distance_3p[1,*] lt distance_3p_req) and (distance_3p[2,*] lt distance_3p_req)
	tell_good_triangle = (triangle_3p[0,*] gt triangle_3p_req) and (triangle_3p[1,*] gt triangle_3p_req) and (triangle_3p[2,*] gt triangle_3p_req)
	tell_good_jvx = (angle_3p_jvx le angle_3p_jvb_req) or (angle_3p_jvx ge 180.-angle_3p_jvb_req)
	tell_good_jvb = (angle_3p_jvb le angle_3p_jvb_req) or (angle_3p_jvb ge 180.-angle_3p_jvb_req)
	tell_noreverse = sign(j3p_regb)*sign(b_ave[0,*]) eq sign(j3p[0,*])
	case dir_suf of
	'': tell_good_dir = tell_good_jvx
	'_regb': tell_good_dir = tell_good_jvb
	'_regbstrict': tell_good_dir = tell_good_jvb and tell_noreverse
	endcase
	tell_good_3p = tell_good_sep3 and tell_good_triangle and tell_good_dir
endif

;;;;;;;; manage data
i_negative_l = where(l_harris lt 0, n_negative_l)
if n_negative_l gt 0 then begin
	l_harris[i_negative_l] = !values.f_nan
	zns[i_negative_l] = !values.f_nan
endif

if ~strcmp(method_suf, '_3p') then begin
	;; optional: mark negative jy to be NaN
	i_negative_jy = where(jy lt 0, n_negative)
	if n_negative gt 0 then begin
		jy[i_negative_jy] = !values.f_nan
	endif
	
	;; compute total current
	I_ttl = 2*Blobe/mu0*1e-9*RE*1000. ;; in A/RE
endif

;;; make the R1-R2 sense current
;; R1 always positive, R2 always negative. jr_regx and jr_regb are for diagnose purpose only.
jr = fltarr(n_elements(jfac))+!values.f_nan
jr_regx = fltarr(n_elements(j3p_regx))+!values.f_nan
jr_regb = fltarr(n_elements(j3p_regb))+!values.f_nan
i_dawn = where(tell_dawn, n_dawn)
if n_dawn gt 0 then begin
	jr[i_dawn] = jfac[i_dawn]
	jr_regx[i_dawn] = j3p_regx[i_dawn]
	jr_regb[i_dawn] = j3p_regb[i_dawn]
endif
i_dusk = where(tell_dusk, n_dusk)
if n_dusk gt 0 then begin
	jr[i_dusk] = -jfac[i_dusk]
	jr_regx[i_dusk] = -j3p_regx[i_dusk]
	jr_regb[i_dusk] = -j3p_regb[i_dusk]
endif

;;;;;;;; plot data
store_data, 'all', data = tell_all
store_data, 'allinx', data = tell_inx
store_data, 'al_low', data = tell_al_low
store_data, 'al_high', data = tell_al_high
if strcmp(method_suf, '_3p') then begin
	store_data, 'good3p', data = tell_good_3p and (tell_dawn or tell_dusk)
	store_data, 'good3p_al_low', data = tell_good_3p and tell_al_low and (tell_dawn or tell_dusk)
	store_data, 'good3p_al_low_limx', data = tell_good_3p and tell_al_low and tell_inx and (tell_dawn or tell_dusk)
	store_data, 'good3p_al_low_dawn', data = tell_good_3p and tell_al_low and tell_dawn
	store_data, 'good3p_al_low_dusk', data = tell_good_3p and tell_al_low and tell_dusk
	store_data, 'good3p_al_low_imfsouth', data = tell_good_3p and tell_al_low and tell_imfsouth
	store_data, 'good3p_al_low_imfnorth', data = tell_good_3p and tell_al_low and tell_imfnorth
	store_data, 'good3p_al_low_dusk', data = tell_good_3p and tell_al_low and tell_dusk
	store_data, 'good3p_al_low_dusk_south', data = tell_good_3p and tell_al_low and tell_dusk and tell_south
	store_data, 'good3p_al_low_dusk_north', data = tell_good_3p and tell_al_low and tell_dusk and tell_north
	store_data, 'good3p_al_high', data = tell_good_3p and tell_al_high and (tell_dawn or tell_dusk)
	store_data, 'good3p_al_high_dawn', data = tell_good_3p and tell_al_high and tell_dawn
	store_data, 'good3p_al_high_dusk', data = tell_good_3p and tell_al_high and tell_dusk
endif

;;;;;;;; choose group of events to plot ;;;;;;;;
;vars_plot = 'allinx' ;; used to reproduce the previous color plots
;vars_plot = 'all'
;vars_plot = 'al_low'
;vars_plot = 'good3p'
;vars_plot = 'good3p_al_low'
;vars_plot = 'good3p_al_low_dawn'
;vars_plot = 'good3p_al_low_dusk'
;vars_plot = 'good3p_al_low_imfsouth'
;vars_plot = 'good3p_al_low_imfnorth'
;vars_plot = 'good3p_al_low_limx'
;vars_plot = 'good3p_al_high'
;vars_plot = 'good3p_al_high_dawn'
;vars_plot = 'good3p_al_high_dusk'
vars_plot = 'good3p_'+['al_low_dawn', 'al_low_dusk'] ;; for paper figure 2
;if strcmp(jvbangle_suf, '_15') then begin
;	vars_plot = 'good3p_'+['al_low', 'al_low'] ;; for double beta plot, supporting matirial
;endif else begin
;	vars_plot = 'good3p_'+['al_low_dawn', 'al_low', 'al_low_dusk', 'al_low', 'al_low'] ;; for paper figure 2
;endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Make location plots based on choice ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;n_minutes = 3 ;; count n minutes difference as different events
;for i = 0, n_elements(vars_plot)-1 do begin
;	get_data, vars_plot[i], data = tell_plot
;	i_plot = where(tell_plot, n_plot)
;	if n_plot gt 0 then begin
;		t_plot = t_all[i_plot]
;		pos_plot = pos[*, i_plot]
;		probes_plot = probes_dual[*, i_plot]
;		rho_plot = rho[i_plot]
;		days_all = strmid(time_string(t_plot), 0, 10)
;		days = days_all[uniq(days_all, sort(days_all))]
;		popen, pic_folder+'/positions'+dir_suf+jvbangle_suf
;		print_options,xsize=3,ysize=8
;		plot, pos_plot[0,*], pos_plot[1,*], /nodata, xtitle = 'X [R!dE!n]', xrange = [-6, -12], ytitle = 'Y [R!dE!n]', title = 'Average Locations', /isotropic
;		;;; find out the pos of each break point of time
;		t_dif = t_plot[1:*]-t_plot[0:-2]
;		i_duan = where((t_dif gt n_minutes*50.) or (t_dif lt 0.), n_dif)
;		if n_dif gt 0 then begin
;			i_duan = [-1, i_duan, -1]
;			for j_duan = 1, n_elements(i_duan)-1 do begin
;				t_this = t_plot[i_duan[j_duan-1]+1:i_duan[j_duan]]
;				pos_this = pos_plot[*, i_duan[j_duan-1]+1:i_duan[j_duan]]
;				probes_this = probes_plot[*, i_duan[j_duan-1]+1:i_duan[j_duan]]
;				oplot, pos_this[0,*], pos_this[1,*]
;				;; diagnose
;				if n_elements(pos_this[0,*]) gt 1 then begin
;					if ~monotonic(pos_this[1,*]) then stop
;				endif
;			endfor
;		endif
;		pclose
;		;;; out put numbers
;		x_plot = pos_plot[0,*]
;		n_use = where(x_plot le -8., n_x_tail)
;		n_use = where(rho_plot ge 8., n_rho_tail)
;		print, 'Number of points:', string(n_plot)
;		print, 'Number of days:', string(n_elements(days))
;		print, 'X<-8 Percentage:', string(n_x_tail/float(n_plot))
;		print, 'Rho>8 Percentage:', string(n_rho_tail/float(n_plot))
;	endif
;endfor
;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Make statistical plots based on choice ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;; choose the regarding quantity (1d plot only)
;qtts = 'habs'
;qtts = ['1obetaFar', '1obetaNear'] ;; for supporting material
qtts = ['habs', 'habs'] ;; for slides
;if strcmp(jvbangle_suf, '_15') then begin
;	qtts = ['habs', '1obeta'] ;; for supporting material, more strict plot
;endif else begin
;	qtts = ['habs', 'znsabs', 'habs', 'habs', '1obeta'] ;; for paper figure 2
;	;qtts = ['1obeta', 'znsabs', '1obeta', '1obeta', '1obeta'] ;; for paper figure 2
;endelse

;;;;;;;;; choose the quantity to check (1d plot only)
qtt2s = replicate('jfac', n_elements(vars_plot)) ;; for proposal
;qtt2s = replicate('jr', n_elements(vars_plot)) ;; for paper

;; design the panel layouts (other parts in folders.pro)
n_panels = n_elements(vars_plot)
if (n_panels ne n_elements(qtts)) or (n_panels ne n_elements(qtts)) then message, 'Vars and qtts not match!'
case n_panels of 
1: begin
	n_p_horiz = 1
	n_p_vert = 1
	size_x = 4.
	size_y = 3.5
	abc = ''
	title_pubs = 'Dawn+Dusk'
	end
;2: begin ;; for paper's supporting material
;	n_p_horiz = 1
;	n_p_vert = 2
;	size_x = 3.5
;	size_y = 6
;	abc = ['(a)', '(b)']
;	title_pubs = ['Dawn+Dusk', '']
;	end
2: begin ;; for slides
	n_p_horiz = 2
	n_p_vert = 1
	space_horiz = 0.005
	size_x = 6
	size_y = 2.5
	abc = ['', '']
	title_pubs = ['Dawn', 'Dusk']
	end
4: begin
	n_p_horiz = 2
	n_p_vert = 2
	size_x = 5
	size_y = 5
	abc = ['(a)', '(b)', '(c)', '(d)']
	end
5: begin
	n_p_horiz = 3
	n_p_vert = 2
	space_horiz = 0.005
	space_vert = 0.135
	if strcmp(plot_type, 'publication') then begin
		size_x = 10.5
		size_y = 7
	endif else begin
		size_x = 8.2
		size_y = 5
	endelse
	abc = ['(a)', '(d)', '(b)', '(c)', '(e)']
	title_pubs = ['Dawn', 'Dawn+Dusk', 'Dusk', 'Dawn+Dusk', 'Dawn+Dusk']
	end
endcase
positions = panel_positions([n_p_horiz, n_p_vert], lr_margins = [left_margin, right_margin], bt_margins = [bot_margin, top_margin], space = [space_horiz, space_vert])

if n_panels eq 5 then begin
	positions[0,1] = positions[0,1]+0.07
	positions[2,1] = positions[2,1]+0.07
	positions[0,5] = positions[0,5]-0.07
	positions[2,5] = positions[2,5]-0.07
	positions = [[positions[*,0:2]], [positions[*,4:5]]]
endif

if strcmp_or(plot_type, ['publication', 'proposal']) then begin
	popen, pic_folder+'/current_stat'+dir_suf+jvbangle_suf
	print_options,xsize=size_x,ysize=size_y
endif

for i = 0, n_elements(vars_plot)-1 do begin
	get_data, vars_plot[i], data = tell_plot
	i_plot = where(tell_plot, n_plot)
	if n_plot gt 0 then begin
		t_plot = t_all[i_plot]
		days_all = strmid(time_string(t_plot), 0, 10)
		days = days_all[uniq(days_all, sort(days_all))]
		ie_plot = i_all[i_plot]
		pos_plot = pos[*, i_plot]
		jfac_plot = jfac[i_plot]
		jr_plot = jr[i_plot]
		beta_plot = beta_ave[i_plot]
		beta_far_plot = beta_far[i_plot]
		beta_near_plot = beta_near[i_plot]
		h_plot = h[i_plot]
		h_far_plot = h_far[i_plot]
		h_near_plot = h_near[i_plot]
		hxy_plot = hxy[i_plot]
		hxy_far_plot = hxy_far[i_plot]
		hxy_near_plot = hxy_near[i_plot]
		l_harris_plot = l_harris[i_plot]
		zns_plot = zns[i_plot] ;; choice 1: first get Z, then averaged. Done in main_current_stat.
		zns_far_plot = zns_far[i_plot] ;; choice 1: first get Z, then averaged. Done in main_current_stat.
		zns_near_plot = zns_near[i_plot] ;; choice 1: first get Z, then averaged. Done in main_current_stat.
	;	zns_plot = l_harris_plot*atanh(h_plot) ;; choice 2: compute Z from average h
		if ~strcmp(method_suf, '_3p') then begin
			jy_plot = jy[i_plot]
			jx_plot = jx[i_plot]
			I_ttl_plot = I_ttl[i_plot]
		endif else begin
			angle_vec_plot = angle_vec[i_plot]
			j3p_regx_plot = j3p_regx[i_plot]
			j3p_regb_plot = j3p_regb[i_plot]
			z_sigma_plot = z_sigma[i_plot]
			l_sigma_plot = l_sigma[i_plot]
		endelse

		;;; for test
	;	stop

;		;;;;;;;;;;;;;;;;;;; make 2d distribution distributions ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;		;;;;;;; horizontal axis
;		;;; Y
;		qtt_horiz_name = 'y'
;		qtt_horiz = transpose(pos_plot[1,*])
;		range_horiz = yrange
;		range_horiz_plot = reverse(yrange)
;		binsize_horiz = binsizey
;		title_horiz = 'Y [RE]'
;
;		;;;;;;; vertical axis
;		;;;; Z
;		;qtt_vert_name = 'z'
;		;qtt_vert = transpose(pos_plot[2,*])
;		;range_vert = zrange
;		;range_vert_plot = zrange
;		;binsize_vert = binsizez
;		;title_vert = 'Z [RE]'
;
;		;;;; h
;		;qtt_vert_name = 'h'
;		;qtt_vert = h_plot
;		;range_vert = [-1., 1.]
;		;range_vert_plot = range_vert
;		;binsize_vert = 0.2
;		;title_vert = 'Bx/Blobe'
;
;		;;; |h|
;		qtt_vert_name = 'habs'
;		qtt_vert = abs(h_plot)
;		range_vert = [0., 1.]
;		range_vert_plot = range_vert
;		binsize_vert = 0.5
;		title_vert = '|Bx/Blobe|'
;
;		;;;;;;; qtt to bin
;		;;; Jx (fac): <0: towards earth, >0: away from earth
;		qtt2_name = 'jx'
;		qtt2bin = jx_plot
;		qtt2_title = 'j!dx!n [nA/m!u2!n]'
;		qtt2_range = jx_range
;
;		;;;; Jy (ctcs)
;		;qtt2_name = 'jy'
;		;qtt2bin = jy_plot
;		;qtt2_title = 'j!dy!n [nA/m!u2!n]'
;		;qtt2_range = jy_range
;
;		;; yz plot
;		bin2d, qtt_horiz, qtt_vert, qtt2bin, xrange = range_horiz, yrange = range_vert, binsize = [binsize_horiz, binsize_vert], binhistogram = counts_2d, averages = qtt2_ave, medians = qtt2_med, stdevs = qtt2_std, xcenters = centers_horiz, ycenters = centers_vert
;
;		;; find out number of events inside each bin and set too few bins to be NaN.
;		n_einbin = events_in_bin(ie_plot, x_arr = qtt_horiz, y_arr = qtt_vert, xcntrs = centers_horiz, ycntrs = centers_vert, counts = counts_2d_ne)
;		i_few = where((counts_2d lt c_pts) or (n_einbin lt c_events), n_few)
;		if n_few gt 0 then begin
;			qtt2_ave[i_few] = !values.f_nan
;			qtt2_med[i_few] = !values.f_nan
;			qtt2_std[i_few] = !values.f_nan
;		endif
;		;;; save data for further use
;		store_data, qtt2_name+'dist'+qtt_horiz_name+qtt_vert_name+'_'+vars_plot[i]+list_suf, data = {type:'2d plot', horiz:centers_horiz, vert:centers_vert, values:qtt2_med, title:vars_plot[i], title_horiz:title_horiz, title_vert:title_vert, title_value:qtt2_title, range_horiz:range_horiz_plot, range_vert:range_vert_plot, range_value:qtt2_range}
;		tplot_save, qtt2_name+'dist'+qtt_horiz_name+qtt_vert_name+'_'+vars_plot[i]+list_suf, filename = save_folder+'/'+qtt2_name+'dist'+qtt_horiz_name+qtt_vert_name+'_'+vars_plot[i]+list_suf
;		;;; plot data
;		plotxyz, centers_horiz, centers_vert, qtt2_med, title=vars_plot[i], xtitle = title_horiz, ytitle = title_vert, ztitle = qtt2_title, /noisotropic, xrange=range_horiz_plot, yrange = range_vert_plot, zrange = qtt2_range
;		makepng, pic_folder+'/'+qtt2_name+'dist'+qtt_horiz_name+qtt_vert_name+'_'+vars_plot[i]+list_suf
;		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		;;;;;;;;;;;;;;;; make 1-D plots ;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;;;; binned quantity
		bin_boundaries = 0
		log1 = 0 ;; default is linear scale
		log2 = 0
		case qtts[i] of
			'Y': begin
				qtt_this = transpose(pos_plot[1,*])
				qtt_range = yrange
				qtt_range_plot = reverse(yrange)
				qtt_title_show = 'Y [RE]'
				binsize_this = binsizey
				vertical = 0
			end
			;;; |h|
			'habs': begin
				qtt_this = abs(h_plot)
				qtt_range = [0., 1.]
				qtt_range_plot = qtt_range
				binsize_this = 0.2 ;; for regB
			;	binsize_this = 0.18 ;; for regBstrict
				qtt_title_show = '|'+jiao_l+'h'+jiao_r+'|'
				vertical = 1
			end
			;;; |hNear|
			'hNearabs': begin
				qtt_this = abs(h_near_plot)
				qtt_range = [0., 1.]
				qtt_range_plot = qtt_range
				binsize_this = 0.2 ;; for regB
			;	binsize_this = 0.18 ;; for regBstrict
				qtt_title_show = '|h!dnear!n|'
				vertical = 1
			end
			;;; |hFar|
			'hFarabs': begin
				qtt_this = abs(h_far_plot)
				qtt_range = [0., 1.]
				qtt_range_plot = qtt_range
				binsize_this = 0.2 ;; for regB
			;	binsize_this = 0.18 ;; for regBstrict
				qtt_title_show = '|h!dfar!n|'
				vertical = 1
			end
			;;; |zns|
			'znsabs': begin
				qtt_this = abs(zns_plot)
				qtt_range = [0., 7.5]
				qtt_range_plot = qtt_range
				binsize_this = 1.5 ;; for regB and regBstrict
				qtt_title_show = '|'+jiao_l+'Z'+jiao_r+minus_sign+'Z!d0!n| [R!dE!n]'
				vertical = 1
				if keyword_set(z_sigma_plot) then begin
					;;; error of L
					rstat, l_sigma_plot, med, lowq, highq
					print, 'Median sigma L:'
					print, [lowq, med, highq]
					;;; error of Z
					rstat, z_sigma_plot, med, lowq, highq
					print, 'Median sigma z:'
					print, [lowq, med, highq]
					stop
				endif
			end
			;;; |znsNear|
			'znsNearabs': begin
				qtt_this = abs(zns_near_plot)
				qtt_range = [0., 10]
				qtt_range_plot = qtt_range
				;binsize_this = 1.5 ;; for regB and regBstrict
				binsize_this = 1. ;; 
				qtt_title_show = '|Z!dnear!n'+minus_sign+'Z!d0!n| [R!dE!n]'
				vertical = 1
			end
			;;; |znsFar|
			'znsFarabs': begin
				qtt_this = abs(zns_far_plot)
				qtt_range = [0., 10]
				qtt_range_plot = qtt_range
				;binsize_this = 1.5 ;; for regB and regBstrict
				binsize_this = 1.
				qtt_title_show = '|Z!dfar!n'+minus_sign+'Z!d0!n| [R!dE!n]'
				vertical = 1
			end
			;;; 1/beta
			'1obeta': begin
				qtt_this = 1/beta_plot
				qtt_range = [0.011, 8]
				if strcmp(jvbangle_suf, '_15') then begin
					n_bin = 6
				endif else begin
					n_bin = 8
				endelse
				rate = (qtt_range[1]/qtt_range[0])^(1./n_bin)
				bin_boundaries = qtt_range[0]*rate^findgen(n_bin+1)
				log1 = 1
				qtt_range_plot = qtt_range
				qtt_title_show = '1/'+jiao_l+beta_letter+jiao_r
				vertical = 1
			end
			'1obetaNear': begin
				qtt_this = 1/beta_near_plot
				qtt_range = [0.011, 8]
				n_bin = 8
				rate = (qtt_range[1]/qtt_range[0])^(1./n_bin)
				bin_boundaries = qtt_range[0]*rate^findgen(n_bin+1)
				log1 = 1
				qtt_range_plot = qtt_range
				qtt_title_show = '1/'+beta_letter+'!dnear!n'
				vertical = 1
			end
			'1obetaFar': begin
				qtt_this = 1/beta_far_plot
				;qtt_range = [0., 10.]
				;binsize_this = 1 ;; for regB

				qtt_range = [0.09, 20]
				n_bin = 8
				rate = (qtt_range[1]/qtt_range[0])^(1./n_bin)
				bin_boundaries = qtt_range[0]*rate^findgen(n_bin+1)
				log1 = 1

				qtt_range_plot = qtt_range
				qtt_title_show = '1/'+beta_letter+'!dfar!n'
				vertical = 1
			end
		endcase

		;;;;;; value quantity
		case qtt2s[i] of
			'l_harris': begin
				qtt2_this = l_harris_plot
				qtt_2_range = [0., 5.]
				qtt_2_title_show = 'l_harris [RE]'
			end
			'i_total': begin
				qtt2_this = I_ttl_plot
				qtt_2_range = [1e5, 6e5]
				qtt_2_title_show = 'i total [A/RE]'
			end
			'jr': begin
				qtt2_this = jr_plot
				;qtt_2_range = [-0.5, 0.5] ;; use this for 2 probes (method_suf = '')
				qtt_2_range = [-0.9, 0.9] ;; use this for 3 probes (method_suf = '_3p'), quiet time
				qtt_2_title_show = 'j!s!d//!r!uR!n (R1'+plus_sign+', R2'+minus_sign+') [nA/m!u2!n]'
			end
			'jfac': begin
				qtt2_this = jfac_plot
				;qtt_2_range = [-0.5, 0.5] ;; use this for 2 probes (method_suf = '')
				qtt_2_range = [-0.9, 0.9] ;; use this for 3 probes (method_suf = '_3p'), quiet time
				qtt_2_title_show = 'j!d//!n [nA/m!u2!n]'
			end
			'angle': begin ;;;; only work for 3 probes (method_suf = '_3p')
				qtt2_this = angle_vec_plot
				i_big = where(qtt2_this gt 90., n_big)
				if n_big gt 0 then qtt2_this[i_big] = 180.-qtt2_this[i_big]
				qtt_2_range = [0., 40.] ;; use this for 3 probes (method_suf = '_3p')
				qtt_2_title_show = 'Angle to the vector [degs]'
			end
		endcase

		;;;;; titles and ticks

		qtt_title_plot = qtt_title_show
		qtt_2_title_plot = qtt_2_title_show
		qtt_ticknames = ''
		qtt_2_ticknames = ''
		if strcmp_or(plot_type, ['publication', 'proposal']) then begin
			noerase = 1
			no_write_pm = 1
			position_this = positions[*,i]
			title = title_pubs[i]
			case n_panels of
				2: begin
					;;;; verticle plot
					;if i eq 0 then begin
					;	if vertical then begin
					;		qtt_2_title_plot = ''
					;		qtt_2_ticknames = replicate(' ', 50)
					;	endif
					;endif else begin
					;	if ~vertical then begin
					;		qtt_2_title_plot = ''
					;		qtt_2_ticknames = replicate(' ', 50)
					;	endif
					;endelse
					;;; horizontal plot
					if i ne 0 then begin
						qtt_title_plot = ''
						qtt_ticknames = replicate(' ', 50)
					endif
				end
				5: begin
					if (i eq 2) or (i eq 3) then begin
						qtt_title_plot = ''
						qtt_ticknames = replicate(' ', 50)
					endif
				end
				else: begin
					noerase = 0
				end
			endcase
		endif else begin
			title = qtt2s[i_qtt2]+' vs '+qtts[i_qtt]+' '+vars_plot[i]
			position_this = 0
			no_write_pm = 0
		endelse

		stat_plot, qtt_this, qtt2_this, k_c = c_pts, bin_range = qtt_range, binsize = binsize_this, qtt_2_range = qtt_2_range, qtt_range = qtt_range_plot, qtt_2_title = qtt_2_title_plot, qtt_title = qtt_title_plot, kinbin = kinbin_all, bincntrs_out = bincntrs, vertical=vertical, qtt_tickname = qtt_ticknames, qtt_2_tickname = qtt_2_ticknames, bin_boundaries = bin_boundaries, title = title, avrg = avrg, std = std, med = med, /no_mean, color_med = 0, color_quar = 0, type_med = 'square', type_quar = 'ebar', position = position_this, noerase = noerase, no_write_pm = no_write_pm, log1 = log1, log2 = log2
		if strcmp_or(qtt2s[i], ['jr', 'jfac']) then begin
			if log1 then begin
				oplot, [0,0.], 10^!y.crange, line = 2
			endif else begin
				oplot, [0,0.], !y.crange, line = 2
			endelse
		endif

		;;; print number of minutes and days
		print, 'Data points:'
		print, n_elements(t_plot)
		print, 'Days involved:'
		print, n_elements(days)

		;;; print figure for normal study
		if strcmp(plot_type, 'test') then begin
			savetitle = qtt2s[i_qtt2]+'_vs_'+qtts[i_qtt]
			makepng, pic_folder+'/'+savetitle+'_'+vars_plot[i]+list_suf+method_suf+dir_suf+jvbangle_suf
		endif else begin
			;;; write abc
			if log1 then begin
				xyouts, !x.crange[0]+0.06*(!x.crange[1]-!x.crange[0]), 10.^(!y.crange[0]+0.9*(!y.crange[1]-!y.crange[0])), abc[i], /data
			endif else begin
				xyouts, !x.crange[0]+0.06*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.9*(!y.crange[1]-!y.crange[0]), abc[i], /data
			endelse
		endelse

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	endif ;; if of n_plot gt 0, have events to plot
endfor
if strcmp_or(plot_type, ['publication', 'proposal']) then pclose

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;; Do statistical study based on events ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; first use only events specified by the variable group
;;; can do only one group of data
;get_data, vars_plot[0], data = tell_plot
;i_plot = where(tell_plot, n_plot)
;
;if n_plot gt 0 then begin
;	t_all = t_all[i_plot]
;	i_all = i_all[i_plot] 
;	pos = pos[*, i_plot]
;	beta_ave = beta_ave[i_plot] 
;	b_ave = b_ave[*, i_plot] 
;	h = h[i_plot]
;	al = al[i_plot]
;	jfac = jfac[i_plot]
;	jr = jr[i_plot]
;	probes_dual = probes_dual[*, i_plot] 
;	if ~strcmp(method_suf, '_3p') then begin
;		;;;; two probes only
;		dxyz = dxyz[*, i_plot]
;		jx = jx[i_plot]
;		jy = jy[i_plot]
;		l_harris = l_harris[i_plot]
;		Blobe = Blobe[i_plot]
;	endif else begin
;		;;; three probes only
;		j3p_regx = j3p_regx[i_plot]
;		j3p_regb = j3p_regb[i_plot] 
;		distance_3p = distance_3p[*, i_plot] 
;		triangle_3p = triangle_3p[*, i_plot] 
;		angle_3p_jvx = angle_3p_jvx[i_plot] 
;		angle_3p_jvb = angle_3p_jvb[i_plot] 
;	endelse
;
;
;	i_uniq = i_all[uniq(i_all, sort(i_all))]
;	
;	h_pre = fltarr(n_elements(i_uniq))+!values.f_nan
;	jr_grwexp = fltarr(n_elements(i_uniq))+!values.f_nan
;	
;	for i = 0, n_elements(i_uniq)-1 do begin
;		k_this = where(i_all eq i_uniq[i], n_this)
;		if n_this lt 1 then message, 'Impossible!'
;		;;; the values of this event
;		t_this = t_all[k_this]
;		pos_this = pos[*, k_this]
;		h_this = h[k_this]
;		al_this = al[k_this]
;		jfac_this = jfac[k_this]
;		jr_this = jr[k_this]
;		l_harris_this = l_harris[k_this] 
;		zns_this = zns[k_this] 
;		if ~strcmp(method_suf, '_3p') then begin
;			dxyz_this = dxyz[*, k_this]
;			jy_this = jy[k_this]
;			jx_this = jx[k_this]
;			Blobe_this = Blobe[k_this]
;		endif else begin
;			j3p_regx_this = j3p_regx[k_this]
;			j3p_regb_this = j3p_regx[k_this]
;		endelse
;	
;		;;;;;;;;;;;;;;;;; do statistics. No need to comment anything, just make sure names are not repeating.
;	
;		;;;;;;;;;;; the h of r-1 cases vs r-2 cases: must use grwexp20 list
;		if t_this[-1]-t_this[0] ge 19*60. then begin 
;			i_pre = where((t_this ge t_this[0]) and (t_this le t_this[0]+5*60.), n_pre) ;; first 5 minutes of the 20-min growth phase range is treated as previous condition
;			i_grwexp = where((t_this ge t_this[-1]-10*60.) and (t_this le t_this[-1]), n_grwexp) ;; 10 minutes is the typical late grw early exp range.
;			if n_pre*n_grwexp gt 0 then begin
;				h_pre[i] = mean(h_this[i_pre], /nan)
;				jr_grwexp[i] = mean(jr_this[i_grwexp], /nan)
;			endif
;		endif
;	endfor
;	
;	;;;;;;;;;;;;;;;; make 1-D plots ;;;;;;;;;;;;;;;;;;;;;;;;;;
;	;;;;;; binned quantity
;	qtt_name = 'habs'
;	qtt_this = abs(h_pre)
;	qtt_range = [0.2, 1.] ;; use this for 3 probes (method_suf = '_3p')
;	binsize_this = 0.2 ;; use this for 3 probes (method_suf = '_3p')
;	;qtt_range = [0.35, 1.] ;; use this for 2 probes (method_suf = '')
;	;binsize_this = 0.2 ;; use this for 2 probes (method_suf = '')
;	qtt_title_show = '|h|'
;	vertical = 1
;	
;	;;;;;; value quantity
;	
;	qtt2_name = 'jrGrwexp'
;	qtt2_this = jr_grwexp
;	;qtt_2_range = [-0.5, 0.5] ;; use this for 2 probes (method_suf = '')
;	qtt_2_range = [-1.5, 1.5] ;; use this for 3 probes (method_suf = '_3p')
;	qtt_2_title_show = 'j R1+ R2- [mV/m]'
;	
;	;;;;; titles
;	title = qtt2_name+' vs '+qtt_name
;	savetitle = qtt2_name+'_vs_'+qtt_name
;	
;	stat_plot, qtt_this, qtt2_this, k_c = c_pts, bin_range = qtt_range, binsize = binsize_this, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title_show, qtt_title = qtt_title_show, kinbin = kinbin_all, bincntrs_out = bincntrs, vertical=vertical, qtt_tickname = qtt_ticknames, qtt_2_tickname = qtt_2_ticknames, bin_boundaries = bin_boundaries, title = title, avrg = avrg, std = std, med = med, /no_mean, color_med = 0, color_quar = 0, type_med = 'square', type_quar = 'ebar'
;	makepng, pic_folder+'/'+savetitle+method_suf+list_suf+dir_suf+jvbangle_suf
;endif else begin
;	message, 'No point satisfying requirement!'
;endelse
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;; Overlap with DFB list ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;listname = 'dfb_list_lead_tail_all'
;events = load_list(listname+'.txt', folder = list_folder)
;n_event = n_elements(events[0,*])
;
;mins_check = 3
;tranges_DFB = [time_double(events[0,*])-mins_check*60., time_double(events[0,*])+mins_check*60.]
;
;;; get all t numbers
;get_data, vars_plot, data = tell_plot
;i_plot = where(tell_plot, n_plot)
;if n_plot gt 0 then begin
;	t_use = t_all[i_plot]
;	h_use = h[i_plot]
;	jfac_use = jfac[i_plot]
;endif
;
;mins_good = lonarr(n_event)
;i_within_plot = -1
;for i = 0, n_event-1 do begin
;	i_within = where((t_use gt tranges_DFB[0,i]) and (t_use lt tranges_DFB[1,i]), n_within)
;	mins_good[i] = n_within
;	if n_within gt 0 then begin
;		i_within_plot = [i_within_plot, i_within]
;	endif
;endfor
;
;i_within_plot = i_within_plot[1:*]
;i_have = where(mins_good gt 0, n_have)
;print, n_have
;
;plot, jfac_use[i_within_plot], abs(h_use[i_within_plot]), psym = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;; for proposal: Find events with j>10 nA/m2 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;i_g3p = where(tell_good_3p)
;jfac_good3 = jfac(i_g3p)
;i_good3 = i_all(i_g3p)
;t_good3 = t_all(i_g3p)
;k_large = where(abs(jfac_good3) gt 5, n_large)
;if n_large gt 0 then begin
;	;i_large = i_all[k_large]
;	;t_large = t_all[k_large]
;	i_large = i_good3[k_large]
;	t_large = t_good3[k_large]
;endif
;
;tranges_print = split_range(t_large, 60*60.)
;
;print, time_string(tranges_print)

stop
end
