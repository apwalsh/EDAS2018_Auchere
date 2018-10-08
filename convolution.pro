!p.multi = 0

loadct, 39

n = 500

signal = fltarr(n)
signal[20:320] = 0.9

psf = fltarr(n)
width = 50.0
psf = exp(-((findgen(n)-n/2)^2.0)/(2.0*width^2.0))/(width*sqrt(2.0*!pi))

;Computation of the TF of the spectrum and PSF
fft_signal = fft(signal)*n
fft_psf = fft(psf)

;Multiplication of the two TFs and inversion
convolution = fft(fft_signal*fft_psf, /inverse)
convolution = abs(convolution)  		 ;Modulus
convolution = shift(convolution, n/2)	 ;and shift to re-center the result

;Display
window, xs=600, ys=600

plot, signal, thick=2, /xs
oplot, psf/max(psf), line=2
oplot, convolution, color=254, thick=2

end
