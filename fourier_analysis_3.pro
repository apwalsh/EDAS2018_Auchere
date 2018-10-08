pro modela, x, a, f, pder

  if a[2] lt 0 then f = 0*x - 1d10 $
  else begin
	  ax = a[0]*x^a[1]
	  f = ax + a[2]
  endelse

  if n_params() ge 4 then $
    pder = [[ax/a[0]], $
		    [ax*alog(x)], $
		    [replicate(1.0d, n_elements(x))]]

end

file = 'light_curve_6.save'

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
psd = psd[1:n/2]

dt = dates[1] - dates[0]
nu_psdt = (1+findgen(n/2))/(n*dt)

;-- Log -log  detrended time-series

frange = [1/(n*dt), 1/(2*dt)]
yrange=[1e-4, 1e4]

plot_oo, nu_psdt, psd, charsize=1, xtitle='frequency', ytitle='Power (sigma_0^2)',$
	/xs, psym=10, title='Time-series - Fourier power spectrum', $
	xr = frange, yr=yrange

weights = 1/smooth(psd, 5, /edge)^2.0
params = [psd[1]*(nu_psdt[1] - nu_psdt[0]), -2, 0]
model_function = 'modela'
fit = curvefit(nu_psdt, psd, weights, params, FUNCTION_NAME=model_function, /DOUBLE, ITMAX=100, TOL=1d-4)

;call_procedure, model_function, nu_psdt, params, fit
oplot, nu_psdt, fit, color=254
oplot, nu_psdt, fit*tresh, color=254, line=2
oplot, nu_psdt, fit*glb_tresh, color=128, line=2


end