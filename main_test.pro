pro main_test
;; test things

del_data, '*'

sc = 'd'
;load_range = ['11 5 20 5', '11 5 20 7'] ;; slow survey range
load_range = ['11 5 20 10', '11 5 20 11'] ;; fast survey range

fill = 1
thm_load_esansst2, probe = sc, trange = load_range, fill = fill
;time_clip, 'th'+sc+'_beta', load_range(0), load_range(1)
tplot, 'th'+sc+'_beta'
end
