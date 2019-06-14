pro main_plot_cartoon
thm_init

dir       = './'
cd, dir
exe_file  = dir + 'magnetosphere_tailcurrent.exe'
file_para  = dir + 'LINE2PARA.txt'
file_trace = dir + 'LINE2TRACE.txt'
file_result= dir + 'LINE2RESULT.txt'
fig        = dir + 'Magnetosphere_tailcurrent'

file_delete,file_trace, /ALLOW_NONEXISTENT
file_delete,file_result, /ALLOW_NONEXISTENT
file_delete,file_para, /ALLOW_NONEXISTENT


;;;; symbol used for plot
usersym, 1.*[-1,0,1,0], 1.*[0,1,0,-1], /fill, thick = 1;; diamond

;;;;;; set plot range
xrange = [2.1, -29]
zrange = [-8, 8]

;;;; settings, no need to change
;        year   iday   ihour  min  isec
para1 = [ 2009,  79,    11,    44,  0] ;; time of equinox
;       vgsex    vgsey    vgsez pdyn   dst   imf_by     imf_bz
para2 = [-304.0,   0.0,   0.0,   3.0, -20.0,  3.0,     -5.0   ]

;;; tail strech ratio
ratio_tails = [0.35, 1.1, 0.35]

;;; construct pos2trace
;        x          y       z      (tail current strengthen ratio).
;;;; dayside field lines, serves all conditions
x_close_dayside = [2., 3.5, 5.5, 7.5, 9.5]
n_close_dayside = n_elements(x_close_dayside)
trace_root_y = 0.

;;; labels and sizes
abc = ['(a)', '(b)', '(c)']
ytitles = ['Quiet Tail', 'Stretched Tail', 'Relaxing Tail']
size_x = 8
size_y = 12
;;; positions
n_p_horiz = 1
n_p_vert = 3
space_horiz = 0.005
space_vert = 0.02
left_margin = 0.1
right_margin = 0.01
bot_margin = 0.01
top_margin = 0.01
positions = panel_positions([n_p_horiz, n_p_vert], lr_margins = [left_margin, right_margin], bt_margins = [bot_margin, top_margin], space = [space_horiz, space_vert])

;;;;;;;; start plotting
popen, 'schematic'
print_options,xsize=size_x,ysize=size_y

for i = 0, n_elements(ratio_tails)-1 do begin
	ratio_tail = ratio_tails[i]
	case ratio_tail of 
		0.35: begin
			;;;; tail closed field lines
			number_close_tail = 9 ;; for atomatic, there will be more manual ones
			inc_tail = 1.15 ;; the rate of ehancement of spacing
			x_close_tail_start = 2.
			inc_times = indgen(number_close_tail-1)+1
			x_close_tail = x_close_tail_start+[0., x_close_tail_start*(1-inc_tail^inc_times)/(1-inc_tail)]
			x_close_tail = [x_close_tail, 39.5, 50, 68, 140]
			x_close_tail = -1*x_close_tail
			;;;; tail lobes
			x_lobe = -1*[0., 6., 13., replicate(20., 4)]
			z_lobe = [15., 17., 18.5, 19., 17., 15., 13.]
			;;;;;; FAC seperation points
			z_sep = 2.5
			z_end = 6.8
		end
		0.28: begin
			;;;; tail closed field lines
			number_close_tail = 9 ;; for atomatic, there will be more manual ones
			inc_tail = 1.15 ;; the rate of ehancement of spacing
			x_close_tail_start = 2.
			inc_times = indgen(number_close_tail-1)+1
			x_close_tail = x_close_tail_start+[0., x_close_tail_start*(1-inc_tail^inc_times)/(1-inc_tail)]
			x_close_tail = [x_close_tail, 36.5, 47, 60, 100]
			x_close_tail = -1*x_close_tail
			;;;; tail lobes
			x_lobe = -1*[0., 6., 13., replicate(20., 4)]
			z_lobe = [15., 17., 18.5, 19., 17., 15., 13.]
			;;;;;; FAC seperation points
			z_sep = 2.5
			z_end = 6.8
		end
		1.1: begin
			;;;;;;; field line trace points
			number_close_tail = 5 ;; for atomatic, there will be more manual ones
			inc_tail = 1.15 ;; the rate of ehancement of spacing
			x_close_tail_start = 2.
			inc_times = indgen(number_close_tail-1)+1
			x_close_tail = x_close_tail_start+[0., x_close_tail_start*(1-inc_tail^inc_times)/(1-inc_tail)]
			x_close_tail = [x_close_tail, 17, 22.5, 29, 36.5, 47, 60, 80, 160]
			x_close_tail = -1*x_close_tail
			;;;; tail lobes
			x_lobe = -1*[0., 6., 13., replicate(20., 6)]
			z_lobe = [15., 17., 18.5, 19., 17., 15., 13., 11, 9]
			;;;;;; FAC seperation points
			z_sep = 0.95
			z_end = 2.3
		end
	endcase

	if i eq 2 then z_sep = z_sep-0.2
	
	n_close_tail = n_elements(x_close_tail)
	n_lobe = n_elements(x_lobe)
	pos2trace = [ $
		; dayside
		[transpose(x_close_dayside), replicate(trace_root_y, 1, n_close_dayside), replicate(0., 1, n_close_dayside), replicate(ratio_tail, 1, n_close_dayside)], $ 
		; nightside closed field lines
		[transpose(x_close_tail), replicate(trace_root_y, 1, n_close_tail), replicate(0., 1, n_close_tail), replicate(ratio_tail, 1, n_close_tail)], $ 
		;; northern lobe
		[transpose(x_lobe), replicate(trace_root_y, 1, n_lobe), transpose(z_lobe), replicate(ratio_tail, 1, n_lobe)], $
	    ;; southern lobe
		[transpose(x_lobe), replicate(trace_root_y, 1, n_lobe), -transpose(z_lobe), replicate(ratio_tail, 1, n_lobe)]]
	
	;;; add the locations for color areas
	scw = [-10., 0, z_sep-0.7, ratio_tail]
	earthward_end = [-6., 0, 0, ratio_tail]
	sep_line = [-10., 0, z_sep, ratio_tail]
	end_line = [-10., 0, z_end, ratio_tail]
	pos2trace = [[pos2trace], [scw], [earthward_end], [sep_line], [end_line]]
	
	
	;;;; write parameters to files for the exe to load
	close,21
	openw,21,file_para
	printf,21,transpose(para1)
	printf,21,transpose(para2)
	close,21
	close,21
	openw,21,file_trace
	printf,21,pos2trace
	close,21
	
	
	while( file_test(file_trace) eq 0) do wait, 0.5
	spawn,exe_file,/NOSHELL, /HIDE
	if(file_test(file_result) ne 1) then message, "Exe file didn't run!"
	nline = file_lines(file_result)
	data  = dblarr(3,nline)
	close,21
	openr,21,file_result
	readf,21,data
	close,21
	
	idx = where(data[0,*] gt 10000.)
	data[*,idx] = !values.d_nan
	n_d = n_elements(data[0,*])
	
	;;; find the three markerlines
	i_nan = where(finite(data[0,*], /nan), n_nan)
	if n_nan gt 3 then begin
		;;; indices of the three field lines
		x_scw_start = transpose(data[0, i_nan[-4]+1:i_nan[-3]-1])
		z_scw_start = transpose(data[2, i_nan[-4]+1:i_nan[-3]-1])
		x_earth_end = transpose(data[0, i_nan[-3]+1:i_nan[-2]-1])
		z_earth_end = transpose(data[2, i_nan[-3]+1:i_nan[-2]-1])
		x_sep = transpose(data[0, i_nan[-2]+1:i_nan[-1]-1])
		z_sep = transpose(data[2, i_nan[-2]+1:i_nan[-1]-1])
		x_high_end = transpose(data[0, i_nan[-1]+1:n_d-1])
		z_high_end = transpose(data[2, i_nan[-1]+1:n_d-1])
		trim_polygon, x_high_end, z_high_end, xrange = xrange
		;;; cut the three lines away from data
		data = data[*, 0:i_nan[-4]-1]
		;;; the two cresents
		x_r2_cres = [x_earth_end, reverse(x_sep)]
		z_r2_cres = [z_earth_end, reverse(z_sep)]
		x_r1_cres = [x_sep, reverse(x_high_end)]
		z_r1_cres = [z_sep, reverse(z_high_end)]
		x_scw = [x_scw_start, reverse(x_sep)]
		z_scw = [z_scw_start, reverse(z_sep)]
		;;;;; trim lines to form smaller middle poligons
		;; R2 current
		x_r2_mid = x_sep
		z_r2_mid = z_sep
		trim_polygon, x_r2_mid, z_r2_mid, xrange = [-8., xrange[1]]
		trim_polygon, x_r2_mid, z_r2_mid, xrange = [-8., -12.]
		;; R1 current
		x_r1_mid = x_high_end
		z_r1_mid= z_high_end
		trim_polygon, x_r1_mid, z_r1_mid, xrange = [-8., xrange[1]]
		trim_polygon, x_r1_mid, z_r1_mid, xrange = [-8., -12.]
	endif else begin
		message, 'No nan point, Cannot happen!'
	endelse
	
	;;;; create plot with empty data
	plot, data[0,*], data[2,*], /iso, xrange = xrange, xstyle = 5, yrange = zrange, ystyle = 5, clip = [xrange[0], zrange[0], xrange[1], zrange[1]], /nodata, position = positions[*,i], /noerase
	
	;;; first, plot the R1 and R2 shapes with poly fill
	;;;; simple plot
	polyfill, x_r1_cres, z_r1_cres, color = 1
	polyfill, x_r1_mid, z_r1_mid, color = 1
	polyfill, x_r2_cres, z_r2_cres, color = 4
	polyfill, x_r2_mid, z_r2_mid, color = 4
	if i eq 2 then begin
		;;; plot SCW
		polyfill, x_scw, z_scw, color = 35
	endif else begin
		;;; draw satellite
		oplot, [-10., !values.f_nan], [1.5, !values.f_nan], psym = 8, symsize = 1.2
	endelse
	;;; then plot field lines
	oplot, data[0,*],data[2,*]
	;;;; below, ranges for diagnose
	;oplot, !x.crange, [z_sep,z_sep]
	;oplot, !x.crange, [z_end,z_end]
	;oplot, !x.crange, [-10,-10.]
	;oplot, !x.crange, [10,10.]
	;oplot, [-8., -8.], !y.crange
	;oplot, [-12., -12.], !y.crange
	;oplot, [-30., -30.], !y.crange
	;;; markers of trace for diagnose
	;plots,pos2trace[0,*],pos2trace[2,*],psym=1
	
	;;;; plot the earth
	y=(findgen(51)-25)/25
	x1=sqrt(1-y^2)
	x2=-sqrt(1-y^2)
	oplot, x1, y
	oplot, x2, y
	oplot, [0,0],[-1,1]
	polyfill, x2, y, color = 0

	;;; write labels
	xyouts, !x.crange[0]-0.06*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.86*(!y.crange[1]-!y.crange[0]), abc[i], /data, charsize = 1.8, align = 0.2
	xyouts, !x.crange[0]-0.06*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.5*(!y.crange[1]-!y.crange[0]), ytitles[i], orientation = -90, /data, charsize = 1.8, align = 0.5
	if i eq 0 then begin
		xyouts, -9.2, -1.8, 'R2', /data, align = 0.5, charsize = 1.6, color = 0
		xyouts, -10, 3.5, 'R1', /data, align = 0.5, charsize = 1.6, color = 255
	endif
	if i eq 2 then begin
		xyouts, -11.54, 0, 'SCW', /data, align = 0.5, charsize = 1.1, color = 255, orientation = 90
	endif
endfor ;; for of i, different panels
pclose

stop
end
