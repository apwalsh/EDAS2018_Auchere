
loadct, 39

!p.multi = [0, 2, 2]

n = 2^12.0  ;number of data points

window, 0

;--------------- M Gaussian Normal (variance 1 then variance V) white noise time series -----------

M = 2

vari = [1.0, 0.5]

for v=0, 1 do begin

	sigma = sqrt(vari[v])

	plot, [-3, 3]*sigma, [0, 1]/sigma, /nodata, /xs, /ys, $
		xtitle='x', ytitle='Probability', $
		title='Variance = ' + strtrim(string(vari[v], format='(F4.1)'), 2)

	bin = sigma/10  ;10 points per sigma

	t = 0
	for i=0, M-1 do begin
		s = sigma*randomn(seed, n, /double)
		t = t + s*s
		h = histogram(s, min=-3*sigma, max=3*sigma, locations=y, bin=bin)/bin/n
		oplot, y, h, psym=10, color=(i+1)*255./M
	endfor

	sigma_t = stddev(t)
	sigma_chi = sqrt(2*M)

	h = histogram(t, min=0, max = 20*sigma_t/sigma_chi, locations=y, bin=bin)/bin/n

	oplot, y, h, psym=10, color=254

	theo_dist = chi2(y/(sigma_t/sigma_chi), M)/(sigma_t/sigma_chi)  ;= exp(-y) for M=2

	oplot, y, theo_dist, line=2

	plot_io, y, h, psym=10, xr=[0, 20*sigma_t/sigma_chi], yr=[1e-3, 1/sigma], color=255, /nodata, $
			/xs, xtitle='x', ytitle='Probability', title='Histogram'

	oplot, y, h, color=255, psym=10
	oplot, y, theo_dist, color=254

endfor

;--------------- Gaussian white noise -------------------

window, 1

print, '------------------ Gaussian white Noise -------------------------'

sigma = 20.0
s = sigma*randomn(seed, n, /double)
s_vari = stddev(s)^2.0
norm = n/s_vari
ff = sqrt(norm)*fft(s, /double)
rff = real_part(ff)
iff = imaginary(ff)

print, 'Variances (signal, real(ft), imaginary(ft)):', s_vari, stddev(rff)^2, stddev(iff)^2

power = rff^2.0d + iff^2.0d

print, 'Power mean:' + string(average(power))
print, 'Power variance:' + string((moment(power))[1])

nu_psdt = dindgen(n)/n

power_th = replicate(1.0, n/2-1)

plot_oo, power[0:n/2-1], xr=[1, n/2], yr=[1e-2, 1e2], ytitle='Power', /xs, $
	title='Gaussian white noise - Power spectrum', psym=10, $
	xtitle='Frequency'
oplot, power_th, color=128

bin = sigma/100.0
histo = histogram(power[0:n/2-1], locations=power_loc, bin=bin, min=0, max=10)/bin
proba = histo/(n/2.0)
plot_io, power_loc, proba, psym=10, yr=[1/n, 10], title='Gaussian white noise - distribution of power', $
	xtitle='Power', xr=[0, 15], ytitle='Probability'

oplot, power_loc, 2*chi2(2*power_loc, 2), color=254 ;= exp(-x)

;AR(1) process

print, '------------------ AR(1) process -------------------------'

alpha = 0.5d

savg = dblarr(n)

s = dblarr(n)
for i=1l, n-1 do s[i] = s[i-1]*alpha + randomn(seed, /double)
s = s * hanning(n)
s_vari = (moment(s))[1]
norm = n/s_vari
power = norm*(abs(fft(s, /double))^2.0d)[0:n/2-1]

power_th = (1 - alpha^2.0)/(1 + alpha^2.0 - 2.0*alpha*cos(2.0*!dpi*dindgen(n/2)/n))

norm_power = power/power_th

print, 'Normalized power mean:' + string(average(norm_power))
print, 'Normalized power variance:' + string((moment(norm_power))[1])

plot_oo, power[0:n/2-1], xr=[1, n/2], ytitle='Power', /ys, psym=10, /xs, $
	title='AR(1) process - Power spectrum', $
	xtitle='Frequency'
oplot, power_th, color=254

bin = sigma/100
histo = histogram(norm_power, locations=power_loc, bin=bin)/bin
proba = histo/(n/2.0)
plot_io, power_loc, proba, psym=10, yr=[1/n, 10], title='AR(1) process - distribution of normalized power', $
	xtitle='Normalized Power', xr=[0, 15], ytitle='Probability'
x = 15*dindgen(100)/100.0
;oplot, x, exp(-x), color=254
oplot, x, 2*chi2(2*x, 2), color=254

end