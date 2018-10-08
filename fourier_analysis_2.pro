file = 'detrended_signal.save'

path = 'C:\Users\Frédéric Auchère\Desktop\ESAC\'

restore, path + file

;--------- Fourier analysis --------------------

window, 0, xs=500, ys=500

!p.multi = [0, 1, 1]

n = n_elements(light_curve)

siglvl = 0.95d
tresh = -alog(1 - siglvl)
glb_tresh = -alog(1 - siglvl^(2.0d/n))

light_curve = light_curve - mean(light_curve)
light_curve = light_curve*hanning(n)
light_curve = light_curve / stddev(light_curve)

norm = n/stddev(light_curve)^2.0
psd = norm*abs(fft(light_curve, /double))^2.0
psd = psd[0:n/2-1]

dt = dates[1] - dates[0]
nu_psdt = (1+findgen(n/2))/(n*dt)

;-- Log -log  detrended time-series

frange = [1/(n*dt), 1/(2*dt)]
yrange=[1e-4, 1e4]

plot_oo, nu_psdt, psd, charsize=1, xtitle='frequency', ytitle='Power (sigma_0^2)',$
	/xs, psym=10, title='Detrended time-series - Fourier power spectrum', $
	xr = frange, yr=yrange
oplot, [1/(n*dt), 2/dt], tresh*[1, 1], line=2, color=254
oplot, [1/(n*dt), 2/dt], glb_tresh*[1, 1], line=2, color=128


light_curve = light_curve*hanning(n)
light_curve = light_curve - mean(light_curve)
light_curve = light_curve / stddev(light_curve)

norm = n/stddev(light_curve)^2.0
psd = norm*abs(fft(light_curve, /double))^2.0
psd = psd[0:n/2-1]


end