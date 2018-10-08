function make_psf, nlambda, dlambda, width

	lambda = dlambda*findgen(nlambda)
	lambda0 = dlambda*nlambda/2.0

	psf = exp(-((lambda-lambda0)^2.0)/(2.0*width^2.0))/(width*sqrt(2.0*!pi))

	return, psf

end