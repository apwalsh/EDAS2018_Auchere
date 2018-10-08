lambda1 = 120.0 ;nanometers
lambda2 = 123.0 ;nanometers
lambda0 = 121.6 ;nanometers
dlambda = 0.01 	;nanometers
intensity = 1.0 ;A.U.
offset = 1.0    ;A.U.
line_width = 0.1;nanometers
psf_width = 0.1 ;nanometers
noise = 0.0   	;noise standard deviation

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

deconvolution = fft(fft_convolution/fft_psf, /inverse, /double)
deconvolution = abs(deconvolution)
deconvolution = shift(deconvolution, nlambda/2)

window, xs=800, ys=400

!p.multi = [0, 2, 1]

loadct, 39

plot, lambda, spectrum, xtitle='Wavelength (nm)', ytitle='Intensity', /xs
oplot, lambda, convolution, line=2
oplot, lambda, psf, line=1
oplot, lambda, deconvolution, line=3, psym=1, color=253

plot_io, abs(fft_convolution), xtitle = 'Frequency Index', ytitle = 'abs(FFT)'
oplot, abs(fft_psf), line=2

!p.multi = 0

end
