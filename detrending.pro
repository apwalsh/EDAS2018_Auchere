file = 'light_curve_6.save'

path = 'C:\Users\Frédéric Auchère\Desktop\ESAC\'

restore, path + file

;parameters
;1: width=4, xrange = [0.2, 0.21]
;2: width=57 xrange=[1.5,3]
;3: width=30, xrange=[40, 41]
;4: width=25, xrange=[55607.075, 55607.09]
;5: width=4, xrange=[2.5, 2.6]
;6  width = 10, xrange=[10, 20]


width = 10
xrange = [10, 20]

loadct, 39
window, 0, xs=500, ys=800



!p.multi = [0, 1, 4]

;----- Original light curve --------------

os = (light_curve - mean(light_curve))/stddev(light_curve)

plot, dates, os, $
	yr=[-3, 3], /xs, xtitle='time', ytitle='signal (sigma)', charsize=2, $
	title='Original signal'

plot, dates, os, $
	psym=1, /xs, xtitle='time', ytitle='signal (sigma)', charsize=2, $
	xr=xrange,title='Original signal zoom' & oplot, dates, os

;--------- Detrending ----------------

detrended_light_curve = light_curve - smooth(light_curve, width, /edge)


ds = (detrended_light_curve - mean(detrended_light_curve))/stddev(detrended_light_curve)
plot, dates, ds, $
	psym=1, yr=[-3, 3], /xs, xtitle='time', ytitle='signal (sigma)', charsize=2, $
	title='Detrended signal' & oplot, dates, ds

plot, dates, ds, $
	psym=1, /xs, xtitle='time', ytitle='signal (sigma)', charsize=2, $
	xr=xrange, title='Detrended signal zoom' & oplot, dates, ds

;---- Save detrended time series -----------------------

light_curve = detrended_light_curve
save, dates, light_curve, file=path + 'detrended_signal.save

end