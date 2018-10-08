function make_spectrum, lambda1, lambda2, dlambda, lambda0, intensity, width, offset, square=square

	nlambda = (lambda2 - lambda1)/dlambda
	lambda = lambda1 + dlambda*findgen(nlambda)

	if keyword_set(square) then begin

		spectrum = replicate(offset, nlambda)
		good = where((lambda gt (lambda0 - width)) and (lambda lt (lambda0 + width)), count)
		if count gt 0 then spectrum[good] = offset + intensity / (count * dlambda)

	endif else $

		spectrum = offset + intensity*exp(-((lambda-lambda0)^2.0)/(2.0*width^2.0))/(width*sqrt(2.0*!pi))

	return, [[lambda], [spectrum]]

end