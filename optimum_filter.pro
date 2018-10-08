pro fit_psd, x, a, f

	f = (a[0] - a[1]*(x/a[2])^2.0d)>a[3]

end

lambda1 = 120.0 ;nanometers
lambda2 = 123.0 ;nanometers
lambda0 = 121.6 ;nanometers
dlambda = 0.01  ;nanometers
intensity = 1.0 ;A.U.
offset = 1.0    ;A.U.
line_width = 0.1;nanometers
psf_width = 0.1 ;nanometers
noise = 0.0     ;noise standard deviation

spectrum = make_spectrum(lambda1, lambda2, dlambda, lambda0, intensity, line_width, offset)

lambda = spectrum[*, 0]
spectrum = spectrum[*, 1]
nlambda = n_elements(lambda)

psf = make_psf(nlambda, dlambda, psf_width)

pad = replicate(0.0, nlambda/2)
psf = [pad, psf, pad]
pad = replicate(offset, nlambda/2)
spectrum = [pad, spectrum, pad]

nlambda = 2*nlambda
lambda = (lambda1 + lambda2)/2.0 - dlambda*nlambda/2.0 + dlambda*findgen(nlambda)

fft_spectrum = fft(spectrum, /double)*nlambda
fft_psf = fft(psf*dlambda, /double)

convolution = fft(fft_spectrum*fft_psf, /double, /inverse)
convolution = abs(convolution)
convolution = shift(convolution, nlambda/2)

lambda = lambda[nlambda/4:3*nlambda/4-1]
spectrum = spectrum[nlambda/4:3*nlambda/4-1]
psf = psf[nlambda/4:3*nlambda/4-1]
convolution = convolution[nlambda/4:3*nlambda/4-1]

nlambda = n_elements(lambda)
convolution = convolution + noise*randomn(seed, nlambda)

fft_convolution = fft(convolution, /double)/nlambda
fft_psf = fft(psf*dlambda, /double)

psd = (abs(fft_convolution)^2.0)[0:nlambda/2]
freq = findgen(nlambda/2+1)

a = fltarr(4)
a[0] = alog10(max(psd[1:*]))
a[1] = 1.0
dum = where(psd gt (max(psd[1:*])/2), count)
a[2] = count*2
a[3] = alog10(mean(psd[0.9*nlambda/2:*]))

weights = replicate(1.0, nlambda)

fit = curvefit(freq[1:*], alog10(psd[1:*]), weights, a, function_name='fit_psd', /noderivative, /double)

noise_model = replicate(10^a[3], nlambda)

a[3] = -100.0d
fit_psd, freq, a, signal_model
signal_model = 10^signal_model
signal_model[0] = psd[0]
signal_model = [signal_model, reverse(signal_model[1:nlambda/2-1])]

optimal_filter = complex(signal_model / (signal_model + noise_model), replicate(0.0d, nlambda), /double)

deconvolution = fft(optimal_filter*fft_convolution/fft_psf, /inverse, /double)
deconvolution = abs(deconvolution)
deconvolution = shift(deconvolution, nlambda/2)

window, xs=800, ys=400

!p.multi = [0, 2, 1]

loadct, 39

plot, lambda, spectrum, xtitle='Wavelength (nm)', ytitle='Intensity', /xs, line=2
oplot, lambda, convolution
oplot, lambda, psf, line=1
oplot, lambda, deconvolution, psym=1, color=253

plot_io, freq, psd
oplot, freq, signal_model, line=2
oplot, freq, noise_model, line=2

!p.multi = 0

end
