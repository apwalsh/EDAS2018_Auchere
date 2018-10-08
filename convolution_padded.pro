!p.multi = 0

loadct, 39

n = 500

signal = fltarr(n)
signal[20:320] = 0.9

npad = n/2
pad = fltarr(npad)
padded_signal = [pad, signal, pad]

psf = fltarr(n)
width = 50.0
psf = exp(-((findgen(n)-n/2)^2.0)/(2.0*width^2.0))/(width*sqrt(2.0*!pi))
padded_psf = [pad, psf, pad]

;Computation of the TF of the spectrum and PSF
fft_signal = fft(padded_signal)*2*n
fft_psf = fft(padded_psf)

;Multiplication of the two TFs and inversion
convolution = fft(fft_signal*fft_psf, /inverse)
convolution = abs(convolution)  		 ;Modulus
convolution = (shift(convolution, n))[npad:npad + n-1]	 ;and shift to re-center the result

;Display
window, xs=600, ys=600

plot, signal, thick=2, /xs
oplot, psf/max(psf), line=2
oplot, convolution, color=254, thick=2

end
