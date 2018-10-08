function chi2, x, k, pdf=pdf

	if keyword_set(pdf) then $
		return, deriv(x, chisqr_pdf(x, k)) $
	else $
		return, (exp(-x/2.0d)*x^(k/2.0-1.0d))/(gamma(k/2.0d)*2.0^(k/2.0d))

end