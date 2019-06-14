pro main_current_diff_plot
;;; load results from main_ctcs_stat_plot and compute the difference from two results.
;;; make sure the two loaded values have the same x and y axis, no interpolation is done here
thm_init
@folders

;;; load data (can load all of them, no matter)
tplot_restore, filename = save_folder+'/'+'jxdistyh_all_grw.tplot'
tplot_restore, filename = save_folder+'/'+'jxdistyh_all_sub.tplot'
tplot_restore, filename = save_folder+'/'+'jxdistyh_all_exp_csv.tplot'
tplot_restore, filename = save_folder+'/'+'jxdistyh_all_rcv_lbl.tplot'

;;; result wll be suf1-suf2
;suf1 = '_exp_csv'
;suf2 = '_rcv_lbl'
suf1 = '_grwexp'
suf2 = '_nongrwexp'

;get_data, 'jxdistyh_all'+suf1, data = jxdist_1
;get_data, 'jxdistyh_all'+suf2, data = jxdist_2
get_data, 'jxdistyhabs_all'+suf1, data = jxdist_1
get_data, 'jxdistyhabs_all'+suf2, data = jxdist_2
;;; values to plot
jx_diff = jxdist_1.values-jxdist_2.values
centers_horiz = jxdist_1.horiz
centers_vert = jxdist_1.vert
title = jxdist_1.title
title_horiz = jxdist_1.title_horiz
title_vert = jxdist_1.title_vert
title_value = jxdist_1.title_value
range_horiz = jxdist_1.range_horiz
range_vert = jxdist_1.range_vert
range_value = jxdist_1.range_value

;;; overwrite settings
range_value = [-0.12, 0.12] ;;for plot

plotxyz, centers_horiz, centers_vert, jx_diff, title=title, xtitle = title_horiz, ytitle = title_vert, ztitle = title_value, /noisotropic, xrange=range_horiz, yrange = range_vert, zrange = range_value
makepng, pic_folder+'/jxdistyh_diff'+suf1+'_minus'+suf2
end
