pro main_event_plot
;;; load the list of events
thm_init
computer='I:'
;computer='/home/jliu'

@folders

pic_folder = pic_folder+'/ctcs'

;; load the list
events = load_trange_list(save_folder+'/ctcs_events.txt', tempfile = save_folder+'/trange_list_template.sav')
;events = events[*, 5]

probes_more = ['a','b','c','d','e']
probes_less = ['a','d','e']

;;;;; make sure the following are the same as main_ctcs_statsettings
;;; range to search
xrange = [-30.,-6.]
rhomax = 12.

;;; sc requirements to compute j
dz_max = 2. ;; in RE, need smaller than this
dxy_ratio2z = 1. ;; in ratio, need smaller than this

;;; al requirement to plot event. The smallest value during the range must be smaller than this to plot event
al_req = -200. 
min_al_after = 10. ;; in minutes. will check all +-this value

;; plot the events
for i = 0, n_elements(events[0,*])-1 do begin
	trange = time_double(events[*,i])
	trange_load = [trange[0]-3.*60., trange[1]+3.*60.]
	trange_al = [trange_load[0], trange_load[1]+min_al_after*60.]
	;;; load data
	del_data, 'thg_idx_al*'
	load_bin_data, trange = trange_al, probes = probes, datatype = 'pseudo_al', datafolder = al_folder, /tclip
	options, 'thg_idx_al', ytitle = 'pseudo AL', ysubtitle = '[nT]'
	get_data, 'thg_idx_al', t, al
	if min(al, /nan) lt al_req then begin
		if trange[0] gt time_double('2010 1 1') then probes = probes_less else probes = probes_more
		del_data, 'th?_state_pos*'
		load_bin_data, trange = trange_load, probes = probes, datatype = 'pos', datafolder = pos_folder, /tclip
		del_data, 'th?_fgs_gsm*'
		load_bin_data, trange = trange_load, probes = probes, datatype = 'fgs', datafolder = fgs_folder, /tclip
		del_data, 'th?_beta*'
		load_bin_data, trange = trange_load, probes = probes, datatype = 'beta', datafolder = beta_folder, /tclip
		del_data, 'th?_Pall*'
		load_bin_data, trange = trange_load, probes = probes, datatype = 'Pall', datafolder = Pall_folder, /tclip

		;;;; generate the time arrays and interpolate all quantities to this
		n_t_arr = fix((trange_load[1]-trange_load[0])/60.)
		t_all = linspace(trange_load[0], trange_load[1], n_t_arr)
		pos_all = fltarr(n_t_arr, n_elements(probes), 3)+!values.f_nan
		beta_all = fltarr(n_t_arr, n_elements(probes))+!values.f_nan
		Pth_all = fltarr(n_t_arr, n_elements(probes))+!values.f_nan
		b_all = fltarr(n_t_arr, n_elements(probes), 3)+!values.f_nan
		;; true or false arrays
		pos_inrange = bytarr(n_t_arr, n_elements(probes))

		for k = 0, n_elements(probes)-1 do begin
			sc = probes[k]
			if tv_exist('th'+sc+'_state_pos_tclip') then begin
			  	get_data, 'th'+sc+'_state_pos_tclip', t_pos, xyz
				xyz_re = xyz/RE
				for i_comp = 0, 2 do begin
					pos_all[*,k,i_comp] = interpol(xyz_re[*,i_comp], t_pos, t_all)
				endfor
				pos_inrange[*,k] = (pos_all[*,k,0] gt xrange[0]) and (pos_all[*,k,0] lt xrange[1]) and (sqrt(pos_all[*,k,1]^2+pos_all[*,k,2]^2) lt rhomax)
			endif
			if tv_exist('th'+sc+'_fgs_gsm_tclip') then begin
				options, 'th'+sc+'_fgs_gsm_tclip', ytitle = 'B '+sc, ysubtitle = '[nT]', labels = ['B!dx', 'B!dy', 'B!dz'], colors = [2, 4, 6]
			  	get_data, 'th'+sc+'_fgs_gsm_tclip', t_b, bxyz
				for i_comp = 0, 2 do begin
					b_all[*,k,i_comp] = interpol(bxyz[*,i_comp], t_b, t_all)
				endfor
			endif
			if tv_exist('th'+sc+'_beta_tclip') then begin
			  	get_data, 'th'+sc+'_beta_tclip', t_beta, beta_this
				beta_all[*,k] = interpol(beta_this, t_beta, t_all)
			endif
			if tv_exist('th'+sc+'_Pall_tclip') then begin
			  	get_data, 'th'+sc+'_Pall_tclip', t_p, Pall
				Pth_all[*,k] = interpol(Pall[*,1], t_p, t_all)
			endif
		endfor ;; for of k, probes

		;;; compute h
		Pttl_2 = nTesla2_to_nPa*(b_all[*,*,0]^2+b_all[*,*,1]^2)+Pth_all
		Bl = sqrt(Pttl_2/nTesla2_to_nPa)
		h_all = b_all[*,*,0]/Bl ;; scale height

		;;; data for plot
		qtts_plot = ''
		qtts_labl = ''
		b_plot = ''

		;;; find suitable ranges and save the current
		for i_sc = 0, n_elements(probes)-2 do begin
			for j_sc = i_sc+1, n_elements(probes)-1 do begin
				d_pos = reform(pos_all[*,i_sc,*]-pos_all[*,j_sc,*], n_t_arr, 3)
				i_good = where((abs(d_pos[*,2]) lt dz_max) and (sqrt(d_pos[*,0]^2+d_pos[*,1]^2) lt dxy_ratio2z*abs(d_pos[*,2])) and pos_inrange[*,i_sc] and pos_inrange[*,j_sc], n_good)
				if n_good gt 0 then begin
					scs = probes[i_sc]+probes[j_sc]
					;; compute current
					d_pos_good = d_pos[i_good, *]
					db_good = reform(b_all[i_good, i_sc, *]-b_all[i_good, j_sc, *], n_good, 3)
					;; add to all result arrays
					t_good = t_all[i_good]
					jy = db_good[*,0]/(mu0*d_pos_good[*,2]*RE*1000)
					jx = db_good[*,1]/(mu0*d_pos_good[*,2]*RE*1000)
					pos_ave = 0.5*reform(pos_all[i_good, i_sc, *]+pos_all[i_good, j_sc, *], n_good, 3)
					pos_dif = abs(reform(pos_all[i_good, i_sc, *]-pos_all[i_good, j_sc, *], n_good, 3))
					beta_ave = 0.5*(beta_all[i_good,i_sc]+beta_all[i_good,j_sc])
					h_ave = 0.5*(h_all[i_good,i_sc]+h_all[i_good,j_sc])
					;; store data
					store_data, 'th'+scs+'_jxy', data = {x: t_good, y: [[jx],[jy]]}
					options, 'th'+scs+'_jxy', ytitle = 'j '+scs, ysubtitle = '[nA/m!u2!n]', labels = ['jx','jy'], colors = [2,4]
					store_data, 'th'+scs+'_beta', data = {x: t_good, y: beta_ave}
					options, 'th'+scs+'_beta', ytitle = 'beta '+scs
					store_data, 'th'+scs+'_h', data = {x: t_good, y: h_ave}
					options, 'th'+scs+'_h', ytitle = 'h '+scs
					;; average position
					store_data, 'th'+scs+'_pos', data = {x: t_good, y: pos_ave}
					split_vec, 'th'+scs+'_pos'
					options, 'th'+scs+'_pos_x', ytitle = 'X'+scs+' [R!dE!n]'
					options, 'th'+scs+'_pos_y', ytitle = 'Y'+scs+' [R!dE!n]'
					options, 'th'+scs+'_pos_z', ytitle = 'Z'+scs+' [R!dE!n]'
					;; position difference
					store_data, 'th'+scs+'_posdif', data = {x: t_good, y: pos_dif}
					split_vec, 'th'+scs+'_posdif'
					options, 'th'+scs+'_posdif_x', ytitle = '|dX'+scs+'| [R!dE!n]'
					options, 'th'+scs+'_posdif_y', ytitle = '|dY'+scs+'| [R!dE!n]'
					options, 'th'+scs+'_posdif_z', ytitle = '|dZ'+scs+'| [R!dE!n]'
					;; add quantities to qtts
					qtts_plot = [qtts_plot, 'th'+scs+'_jxy', 'th'+scs+'_beta']
					qtts_labl = [qtts_labl, 'th'+scs+'_pos_x', 'th'+scs+'_pos_y', 'th'+scs+'_pos_z', 'th'+scs+'_h', 'th'+scs+'_posdif_x', 'th'+scs+'_posdif_y', 'th'+scs+'_posdif_z']
					;; add quantities to b_plot
					if ~strcmp_or('th'+probes[i_sc]+'_fgs_gsm_tclip', b_plot) then b_plot = [b_plot, 'th'+probes[i_sc]+'_fgs_gsm_tclip']
					if ~strcmp_or('th'+probes[j_sc]+'_fgs_gsm_tclip', b_plot) then b_plot = [b_plot, 'th'+probes[j_sc]+'_fgs_gsm_tclip']
				endif
			endfor ;; for of i_sc
		endfor ;; for of j_sc

		;;; make plot
		if n_elements(b_plot) gt 1 then begin
			b_plot = b_plot[1:*]
			sc_used = strmid(b_plot, 2, 1)
			split_vec, b_plot
			components = ['x', 'y', 'z']
			for i_xyz = 0, n_elements(components)-1 do begin
				store_data, 'all_b'+components[i_xyz], data = 'th?_fgs_gsm_tclip_'+components[i_xyz]
				options, 'all_b'+components[i_xyz], ytitle = 'B!d'+components[i_xyz], ysubtitle = '[nT]'
			endfor
			for i_sc = 0, n_elements(sc_used)-1 do begin
				options, 'th'+sc_used[i_sc]+'_fgs_gsm_tclip_?', colors = thm_probe_color(sc_used[i_sc]), labels = sc_used[i_sc]
			endfor
			ball_plot = ['all_bx', 'all_by', 'all_bz']
		endif else ball_plot = ''
		if n_elements(qtts_plot) gt 1 then qtts_plot = qtts_plot[1:*]
		if n_elements(qtts_labl) gt 1 then qtts_labl = qtts_labl[1:*]
		tplot, ['thg_idx_al_tclip', ball_plot, qtts_plot], var_label = reverse(qtts_labl), trange = trange_load, title = 'event '+strcompress(string(i),/remove)+' '+time_string(mean(trange), format = 6, precision = -1)
		timebar_mass, 0, varname = [ball_plot, qtts_plot], /databar, line = 1
		makepng, pic_folder+'/event'+strcompress(string(i),/remove)+'_'+time_string(mean(trange), format = 6, precision = -1)
	endif ;; if of AL meet requirement
endfor
stop
end
