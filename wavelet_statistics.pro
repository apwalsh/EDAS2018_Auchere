pro chi2_fit, x, a, f

	if a lt 0 then f = replicate(-1d20, n_elements(x)) $
	else f = chi2(x*a, a, /pdf)*a

end

loadct, 39

n = 2^7
dj = 8
ntests = 1000

mother = 'Morlet'

nsigmas = 20
sigma_max = 20.0
sigmas =  sigma_max*findgen(nsigmas)/nsigmas

pad = 1
exclude_coi = 1

dj = 1.0/8
noct = alog(n)/alog(2.0) - 2
j1 = noct/dj
dt = 1.0
s0 = 2.0*dt
k0 = 6.0

;--- Morlet --------------------------
dofmin = 2.0
morlet_coi = 0.730472
gamma_m = 2.32 ;empirical factor for time averaging Morlet
morlet_dj0 = 0.6 ;empirical factor for scale averaging Morlet
fourier_factor = (4*!pi)/(k0 + sqrt(2+k0^2))

wvl_avg = fltarr(ntests)

sst = randomn(seed, n, /double)
wave = WAVELET(sst,dt,PERIOD=period,SCALE=scale,S0=s0, $
                PAD=pad,COI=coi,DJ=dj,J=j1,MOTHER=mother)

mask = replicate(1.0, n, n_elements(scale))

if exclude_coi eq 1 then begin
	for i=0, n_elements(scale)-1 do begin
		mask[0:(scale[i]/morlet_coi/dt)<(n/2), i] = 0.0
		mask[n-((scale[i]/morlet_coi/dt)<(n/2)):*, i] = 0.0
	endfor
endif

pdf = fltarr(nsigmas)
scale_pdf = fltarr(nsigmas)

wvl_spectra = dblarr(ntests, j1+1)

if exclude_coi eq 1 then norm_coi = 1.0 - 2.0*scale/morlet_coi/n $
else norm_coi = replicate(1.0, j1+1)
gd_norm_coi = where(norm_coi gt 0)

if exclude_coi eq 1 then valid_scales = where((scale/morlet_coi) le (n/2.0), nvalid_scales) $
else begin
	valid_scales = findgen(j1+1)
	nvalid_scales = j1+1
endelse

tot_mask = total(mask, /double)

avg_spectrum = fltarr(n, j1+1)

sigmas_scales = fltarr(nsigmas, j1+1)
probas = exp(-sigmas) ;local proabability associated to sigma for a chi2 distribution

na = n - scale/(fourier_factor/sqrt(2))
dof = dofmin*SQRT( 1 + (na*dt/gamma_m/scale)^2.0d ) ; [Eqn(23)]
dof = dof > dofmin

for iscl=0, j1 do $  ;corresponding local probabilities for the time averaged spectrum
	for isig =0, nsigmas-1 do $
		sigmas_scales[isig, iscl] = CHISQR_CVF(probas[isig], dof[iscl])/dof[iscl]

for it=0L, ntests-1 do begin

	sst = randomn(seed, n, /double)

    wave = WAVELET(sst,dt,param=k0,PERIOD=period,SCALE=scale,S0=s0, $
                PAD=pad,COI=coi,DJ=dj,J=j1,MOTHER=mother)

	power = (abs(wave))^2.0d

	avg_spectrum = avg_spectrum + power

	power = power*mask

	dum = total(power, 1, /double)/n
	dum[gd_norm_coi] = dum[gd_norm_coi] / norm_coi[gd_norm_coi]

    wvl_spectra[it, *] = dum

	wvl_avg[it] = total(power, /double)/tot_mask

	maxi = max(power)
	tresh = round(nsigmas*(maxi/sigma_max))<(nsigmas-1)
	pdf[0:tresh] = pdf[0:tresh] + 1.0

	isig = nsigmas
	repeat begin

		isig = isig - 1
		if exclude_coi eq 1 then $
			above = where(reform(wvl_spectra[it, valid_scales]) ge reform(sigmas_scales[isig, valid_scales]), count) $
		else $
        	above = where(reform(wvl_spectra[it, *]) ge reform(sigmas_scales[isig, *]), count)

	endrep until count ge 1

	scale_pdf[0:isig] = scale_pdf[0:isig] + 1

endfor

avg_spectrum = avg_spectrum / ntests

pdf = pdf / ntests
scale_pdf = scale_pdf / ntests

;------------------------------------------------------------------------------------------------

histo_min = 0.0
histo_max = 10.0
histo_bin = 0.01
histo_nbins = (histo_max - histo_min)/histo_bin


window, xs=900, ys=600
!p.multi = [0, 3, 2]

;-------- Distribution of average power in wavelet spectra -----------------------

fdof = n_elements(power)*dj/3.0d
histo = histogram(wvl_avg, min=histo_min, bin=histo_bin, nbins=histo_nbins, $
	locations=histo_x, /nan, /l64)
histo = histo/total(histo)/histo_bin
good = where(histo gt 0)
weights = 1.0/histo[good]
res = curvefit(histo_x[good], histo[good], weights, fdof, /noderivative, $
	function_name='chi2_fit', /double)

plot, histo_x[good], histo[good], psym=10, charsize=2, $
	xtitle='Power (sigma)', ytitle='Counts', title='mean power in wavelet spectra', $
	subtitle='DOF=' + strtrim(string(fdof, format='(F5.1)'), 2)
oplot, histo_x[good], res, color=254, thick=2


;-------- Distribution of power at a given scale in time-averaged wavelet spectra -----------------------

iscale = 30  ;just pick one scale
fdof = dof[iscale]
histo = histogram(wvl_spectra[*, iscale], min=histo_min, bin=histo_bin, $
	nbins=histo_nbins, locations=histo_x, /nan, /l64)
histo = histo/total(histo)/histo_bin
good = where(histo gt 0)
weights = 1.0/histo[good]
res = curvefit(histo_x[good], histo[good], weights, fdof, /noderivative, $
	function_name='chi2_fit', /double)

plot, histo_x[good], histo[good], psym=10, charsize=2, $
	xtitle='Power (sigma)', ytitle='Counts', title='t-averaged power for scale ' + strtrim(string(iscale), 2), $
	subtitle='DOF=' + strtrim(string(fdof, format='(F5.1)'), 2) $
					+ ' (TC98=' + strtrim(string(dof[iscale], format='(F5.1)'), 2) + ')'
oplot, histo_x[good], res, color=254, thick=2

;----------- TC98 fit -----------------------
oplot, histo_x[good], dof[iscale]*chi2(histo_x[good]*dof[iscale], dof[iscale]), color=128, thick=2

;-------- Distribution of average power in time-averaged wavelet spectra -----------------------

fdof = total(dof)*dj
histo = histogram(total(wvl_spectra[*, valid_scales], 2)/nvalid_scales, min=histo_min, $
	bin=histo_bin, nbins=histo_nbins, locations=histo_x, /nan, /l64)
histo = histo/total(histo)/histo_bin
good = where(histo gt 0)
weights = 1.0/histo[good]
res = curvefit(histo_x[good], histo[good], weights, fdof, /noderivative, $
	function_name='chi2_fit', /double)

plot, histo_x[good], histo[good], psym=10, charsize=2, $
	xtitle='Power (sigma)', ytitle='Counts', title='mean power in t-averaged spectrum', $
	subtitle='DOF=' + strtrim(string(fdof, format='(F5.1)'), 2)
oplot, histo_x[good], res, color=254, thick=2

;-------- Global vs. local probabilities in the wavelet sepctra -----------------------

plot_oo, probas, pdf, xtitle='P(m)', ytitle='P_g(m)', charsize=2, xr=[1e-8, 1], yr=[1e-4, 1], $
	title='Wavelet: global vs. local probas'

aexp = 0.810*(nvalid_scales*n*dj)^0.011d			;Eq. 16
nexp = 0.491*(nvalid_scales*n*dj)^0.926d			;Eq. 17

glb_probas = 1 - (1 - probas^aexp)^nexp

oplot, probas, glb_probas, color=254, line=2

;-------- Global vs. local probabilities in the time-averaged wavelet sepctra -----------------------

plot_oo, probas, scale_pdf, xtitle='P(m)', ytitle='P_g(m)', charsize=2, xr=[1e-8, 1], yr=[1e-4, 1], $
	title='t-averaged wavelet: global vs. local probas'

aexp_scl = 0.805d + 0.45*2.0^(-nvalid_scales*dj)	;Eq. 19	;
nexp_scl = 1.136*(nvalid_scales*dj)^1.2d			;Eq. 20

glb_probas = 1 - (1 - probas^aexp_scl)^nexp_scl

oplot, probas, glb_probas, color=254, line=2


end
