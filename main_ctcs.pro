pro main_ctcs
;; compute cross-tail current
;; Note: in the cdfs, for a given geomagnetic level, X,Y,Z are the same.
thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

;; the location to check the cross tail current distribution. A single number for interpolation to a plane, a range for summing all data points falling into the range.
x_check = [-10.] ;; for interpolation
;x_check = [-9.5, -10.5] ;; for summing this range
nypts = 50 ;; number of points in Y to check
nzpts = 30 ;; number of points in Z for interpolation

;; the beta value to mark the boundary of the cross tail current sheet
beta_psbl = 0.5

;; number of points in regular grid
nx_reg = 100
ny_reg = 50
nz_reg = 25

folder_cdf = 'ctcs_values'

act_sufs = ['_pd=2_sunspot=50_100', '_pd=2_sunspot=50_8000']

for i = 0, n_elements(act_sufs)-1 do begin
	filename_c = folder_cdf+'/'+'currents'+act_sufs[i]+'.cdf'
	filename_b = folder_cdf+'/'+'mag_field'+act_sufs[i]+'.cdf'
	filename_p = folder_cdf+'/'+'pressure_isotropic'+act_sufs[i]+'.cdf'
	;; current
	cdf_var_get, filename_c, varnames = ['xFull', 'yFull', 'zFull', 'jCross', 'jPar'], /ncdf, suffix = act_sufs[i]
	;; pressure
	cdf_var_get, filename_p, varnames = 'presEq', /ncdf, suffix = act_sufs[i]
	;; bfield
	cdf_var_get, filename_b, varnames = 'bFull', /ncdf, suffix = act_sufs[i]
	;;; compute beta
	get_data, 'presEq'+act_sufs[i], data = Peq
	Peq = Peq[*,1:*]
	Peq = reform(Peq, 1, n_elements(Peq[*, 0]), n_elements(Peq[0, *]))
	Pth = rebin(Peq, 37, n_elements(Peq[0, *, 0]), n_elements(Peq[0, 0, *]))
	get_data, 'bFull'+act_sufs[i], data = Bmag
	Bmag = Bmag[37:73,*,1:*] ;; 37:73: only half of the points, because north and south are symmetric.
	beta_value = Pth*1e-9/((Bmag*1e-9)^2/(2*mu0))
	store_data, 'beta'+act_sufs[i], data = beta_value
endfor

;;; process the location values
loc_names = ['xFull'+act_sufs[0], 'yFull'+act_sufs[0], 'zFull'+act_sufs[0], 'xFull'+act_sufs[1], 'yFull'+act_sufs[1], 'zFull'+act_sufs[1]]
for i = 0, n_elements(loc_names)-1 do begin
	get_data, loc_names[i], data = loc
	loc = loc[37:73,*,1:*]
	loc_min = min(loc, /nan)
	loc_max = max(loc, /nan)
	loc_diff = loc-loc_min
	loc_normal = loc_diff/(max(loc_diff, /nan)*(1.+1e-5))
	if strmatch(loc_names[i], 'x*') then begin
		arr_regular = linspace(loc_min, loc_max, nx_reg)
		mat_regular = rebin(arr_regular, nx_reg, ny_reg, nz_reg)
	endif
	if strmatch(loc_names[i], 'y*') then begin
		arr_regular = transpose(linspace(loc_min, loc_max, ny_reg))
		mat_regular = rebin(arr_regular, nx_reg, ny_reg, nz_reg)
	endif
	if strmatch(loc_names[i], 'z*') then begin
		arr_regular = linspace(loc_min, loc_max, nz_reg)
		arr_regular = reform(arr_regular, 1, 1, nz_reg)
		mat_regular = rebin(arr_regular, nx_reg, ny_reg, nz_reg)
	endif
	store_data, loc_names[i]+'_processed', data = {data:loc, min:loc_min, max:loc_max, normal:loc_normal, regular_mat:mat_regular, regular_arr:arr_regular, n:n_elements(loc)}
endfor

;;;;; compute the total current
;;;; reduce the data to be half; faster
get_data, 'xFull'+act_sufs[0]+'_processed', data = x_c
get_data, 'xFull'+act_sufs[1]+'_processed', data = x_c_g
;;;; determine Y range
get_data, 'yFull'+act_sufs[0]+'_processed', data = y_c
get_data, 'yFull'+act_sufs[1]+'_processed', data = y_c_g
yrange = [min([y_c.min, y_c_g.min]), max([y_c.max, y_c_g.max])]
;;;; determine Z range
get_data, 'zFull'+act_sufs[0]+'_processed', data = z_c
get_data, 'zFull'+act_sufs[1]+'_processed', data = z_c_g
zrange = [min([z_c.min, z_c_g.min]), max([z_c.max, z_c_g.max])]
;;;; current data
;;; cross tail current
get_data, 'jCross'+act_sufs[0], data = j_q
j_q = -j_q[37:73,*,1:*]
get_data, 'jCross'+act_sufs[1], data = j_g
j_g = -j_g[37:73,*,1:*]
;;; parallel current
get_data, 'jPar'+act_sufs[0], data = fac_q
fac_q = fac_q[*,1:*]
get_data, 'jPar'+act_sufs[1], data = fac_g
fac_g = fac_g[*,1:*]
;; facs are 2-d because they are along with field lines. Rebin the data to easier interpole later.
fac_q = rebin(fac_q, n_elements(fac_q[*,0]), n_elements(fac_q[0,*]), n_elements(j_q[*,0,0]))
fac_q = transpose(fac_q, [2,0,1])
fac_g = rebin(fac_g, n_elements(fac_g[*,0]), n_elements(fac_g[0,*]), n_elements(j_g[*,0,0]))
fac_g = transpose(fac_g, [2,0,1])
;;;; beta value
get_data, 'beta'+act_sufs[0], data = beta_q
get_data, 'beta'+act_sufs[1], data = beta_g

;;; 3-d interpolation
y_check = linspace(yrange[0], yrange[1], nypts)

if n_elements(x_check) lt 2 then begin
	;;; interpolation
	;; first, adjust to even grids
	j_q_regular = interp_cic(j_q, x_c.normal*nx_reg, nx_reg, y_c.normal*ny_reg, ny_reg, z_c.normal*nz_reg, nz_reg, /average)
	j_g_regular = interp_cic(j_g, x_c_g.normal*nx_reg, nx_reg, y_c_g.normal*ny_reg, ny_reg, z_c_g.normal*nz_reg, nz_reg, /average)
	fac_q_regular = interp_cic(fac_q, x_c.normal*nx_reg, nx_reg, y_c.normal*ny_reg, ny_reg, z_c.normal*nz_reg, nz_reg, /average)
	fac_g_regular = interp_cic(fac_g, x_c_g.normal*nx_reg, nx_reg, y_c_g.normal*ny_reg, ny_reg, z_c_g.normal*nz_reg, nz_reg, /average)
	beta_q_regular = interp_cic(beta_q, x_c.normal*nx_reg, nx_reg, y_c.normal*ny_reg, ny_reg, z_c.normal*nz_reg, nz_reg, /average)
	beta_g_regular = interp_cic(beta_g, x_c_g.normal*nx_reg, nx_reg, y_c_g.normal*ny_reg, ny_reg, z_c_g.normal*nz_reg, nz_reg, /average)
	z_q_regular = interp_cic(z_c.data, x_c.normal*nx_reg, nx_reg, y_c.normal*ny_reg, ny_reg, z_c.normal*nz_reg, nz_reg, /average)
	z_g_regular = interp_cic(z_c_g.data, x_c_g.normal*nx_reg, nx_reg, y_c_g.normal*ny_reg, ny_reg, z_c_g.normal*nz_reg, nz_reg, /average)
	;;;; sum over Z to get I distribution as function of X and Z
	;; quiet
	i_outpsbl_q = where(beta_q_regular lt beta_psbl, n_outpsbl_q)
	if n_outpsbl_q gt 0 then begin
		j_q_regular[i_outpsbl_q] = !values.f_nan
		z_q_regular[i_outpsbl_q] = !values.f_nan
	endif
	I_q = total(j_q_regular,3,/nan)
	z_psbl_q = max(z_q_regular, dimension = 3, /nan)
	i_q_tail = where(x_c.regular_arr lt 0, n_q_tail)
	if n_q_tail gt 0 then begin 
		I_q = I_q[i_q_tail, *]
		z_psbl_q = z_psbl_q[i_q_tail, *]
		x_c_q_tail = x_c.regular_arr[i_q_tail]
	endif
	;; growth
	i_outpsbl_g = where(beta_g_regular lt beta_psbl, n_outpsbl_g)
	if n_outpsbl_g gt 0 then begin
		j_g_regular[i_outpsbl_g] = !values.f_nan
		z_g_regular[i_outpsbl_g] = !values.f_nan
	endif
	z_psbl_g = max(z_g_regular, dimension = 3, /nan)
	I_g = total(j_g_regular,3,/nan)
	i_g_tail = where(x_c_g.regular_arr lt 0, n_g_tail)
	if n_g_tail gt 0 then begin
		I_g = I_g[i_g_tail, *]
		z_psbl_g = z_psbl_g[i_g_tail, *]
		x_c_g_tail = x_c_g.regular_arr[i_g_tail]
	endif
	;;; re-mesh the grid
	;x_c_q_tail_mat = rebin(x_c_q_tail, n_q_tail, y_c.n)
	;x_c_g_tail_mat = rebin(x_c_g_tail, n_g_tail, y_c_g.n)
	;y_c_q_tail_mat = rebin(y_c.regular_arr, n_q_tail, y_c.n)
	;y_c_g_tail_mat = rebin(y_c_g.regular_arr, n_g_tail, y_c_g.n)
	;;; then interpolate
	;j_q_result = grid3(x_c.regular_mat, y_c.regular_mat, z_c.regular_mat, j_q_regular, x_check, y_check, /grid, dtol = 0)
	;j_g_result = grid3(x_c_g.regular_mat, y_c_g.regular_mat, z_c_g.regular_mat, j_g_regular, x_check, y_check, /grid, dtol = 0)
	;fac_q_result = grid3(x_c.regular_mat, y_c.regular_mat, z_c.regular_mat, fac_q_regular, x_check, y_check, /grid, dtol = 0)
	;fac_g_result = grid3(x_c_g.regular_mat, y_c_g.regular_mat, z_c_g.regular_mat, fac_g_regular, x_check, y_check, /grid, dtol = 0)

	;;; compute the current
	I_q_mid = fltarr(n_elements(y_c.regular_arr))
	z_psbl_q_mid = fltarr(n_elements(y_c.regular_arr))
	I_g_mid = fltarr(n_elements(y_c_g.regular_arr))
	z_psbl_g_mid = fltarr(n_elements(y_c_g.regular_arr))
	for k = 0, n_elements(y_c.regular_arr)-1 do begin
		I_q_mid[k] = interpol(I_q[*,k], x_c_q_tail, x_check)
		z_psbl_q_mid[k] = interpol(z_psbl_q[*,k], x_c_q_tail, x_check)
	endfor
	for k = 0, n_elements(y_c_g.regular_arr)-1 do begin
		I_g_mid[k] = interpol(I_g[*,k], x_c_g_tail, x_check)
		z_psbl_g_mid[k] = interpol(z_psbl_g[*,k], x_c_g_tail, x_check)
	endfor
	I_quiet = interpol(I_q_mid, y_c.regular_arr, y_check)
	I_growth = interpol(I_g_mid, y_c_g.regular_arr, y_check)
	z_psbl_q = interpol(z_psbl_q_mid, y_c.regular_arr, y_check)
	z_psbl_g = interpol(z_psbl_g_mid, y_c_g.regular_arr, y_check)
endif else begin
	;;; collect values that falls into each Y range
	I_quiet = fltarr(nypts)
	I_growth = fltarr(nypts)
	z_q = z_c.data
	z_g = z_c_g.data
	;; mark all data out of the X range or out of the psbl (boundary of the ctcs) to be NaN
	i_out_q = where((x_c.data gt max(x_check)) or (x_c.data lt min(x_check)) or (beta_q lt beta_psbl), n_out_q)
	i_out_g = where((x_c_g.data gt max(x_check)) or (x_c_g.data lt min(x_check)) or (beta_q lt beta_psbl), n_out_g)
	if n_out_q gt 0 then begin
		j_q[i_out_q] = !values.f_nan
		fac_q[i_out_q] = !values.f_nan
		z_q[i_out_q] = !values.f_nan
	endif
	if n_out_g gt 0 then begin
		j_g[i_out_g] = !values.f_nan
		fac_g[i_out_g] = !values.f_nan
		z_g[i_out_g] = !values.f_nan
	endif
	;; collect the points regarding Y
	for i = 0, nypts-1 do begin
		;; the range for searching values
		if i eq 0 then ydiff = y_check[i+1]-y_check[i] else ydiff = y_check[i]-y_check[i-1]
		yrange_search = [y_check[i]-0.5*ydiff, y_check[i]+0.5*ydiff]
		i_inrange_q = where((y_c.data gt min(yrange_search)) and (y_c.data lt max(yrange_search)), n_inrange_q)
		i_inrange_g = where((y_c_g.data gt min(yrange_search)) and (y_c_g.data lt max(yrange_search)), n_inrange_g)
		;; sum the values
		if n_inrange_q gt 0 then begin
			I_quiet[i] = total(j_q[i_inrange_q], /nan)
			z_psbl_q[i] = max(z_q[i_inrange_q],/nan)
		endif else begin
			I_quiet[i] = !values.f_nan
			z_psbl_q[i] = !values.f_nan
		endelse
		if n_inrange_g gt 0 then begin
			I_growth[i] = total(j_g[i_inrange_g], /nan)
			z_psbl_g[i] = max(z_g[i_inrange_g],/nan)
		endif else begin
			I_growth[i] = !values.f_nan
			z_psbl_g[i] = !values.f_nan
		endelse
	endfor
endelse

dI = I_growth-I_quiet

;;;;;;;;;;;;; plot results
;;; plot the I distribution
;plot, y_check, I_growth, /nodata, xtitle = 'Y [RE]', ytitle = 'I', xrange = [yrange[1], yrange[0]], xstyle = 1, yrange = [0., max(I_growth)*1.1], ystyle = 1, title = 'b-quiet, r-growth'
;oplot, y_check, I_quiet, color = 2
;oplot, y_check, I_growth, color = 6
;
;;; plot the cs thickness
;plot, y_check, z_psbl_q, /nodata, xtitle = 'Y [RE]', ytitle = 'Z boundary [RE]', xrange = [yrange[1], yrange[0]], xstyle = 1, yrange = [0., max(z_psbl_q)*1.1], ystyle = 1, title = 'b-quiet, r-growth'
;oplot, y_check, z_psbl_q, color = 2
;oplot, y_check, z_psbl_g, color = 6

;; plot the distribution of j//
help, fac_q_regular
help, fac_g_regular
help, y_c.regular_arr
help, y_c_g.regular_arr
help, z_c.regular_arr
help, z_c_g.regular_arr
;plotxyz, fac_q_regular
;fac_g_regular
;fac_g_regular-fac_q_regular

stop
end
