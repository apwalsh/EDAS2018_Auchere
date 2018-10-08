;Converts the global confidence level (vertical axes of Figures 6 & 7) to the corresponding local
;confidence level (Horizontal axes of Figures 6 & 7)
;Npoints is the number of data points outside the COI in the wavelet spectrum or in the time-averaged wavelet spectrum
;Specify /gws if you want to compute confidence levels for the time-averaged wavelet spectrum

function global_to_local_siglvl, siglvl, Mother, Npoints, dj, gws=gws

	case Mother of
		'Morlet': begin
			aexp = 0.810*(npoints*dj)^0.011d			;Eq. 16
			nexp = 0.491*(npoints*dj)^0.926d			;Eq. 17

			aexp_scl = 0.805d + 0.45*2.0^(-npoints*dj)	;Eq. 19	;
			nexp_scl = 1.136*(npoints*dj)^1.2d			;Eq. 20
		end
		;!!!!!Warning !!! The coefficients have been derived for the Paul wavelet as for the Morlet wavelet,
		;but I did not test the results as extensively.
		'Paul': begin
			aexp = 0.817*(npoints*dj)^0.011d			;Eq. 21
			nexp = 0.320*(npoints*dj)^0.926d			;Eq. 22

			aexp_scl = 1.02d + 0.70*21.42^(-npoints*dj)	;Eq. 23
			nexp_scl = 1.0 + 0.56*(npoints*dj)^1.2d		;Eq. 24
		end

	endcase

	if ~keyword_set(gws) then $
		return,  1.0d - (1.0d - siglvl^(1.0/nexp))^(1.0/aexp) $ ;Eq. 18 (anydtperiod_siglvl)
	else $
		return, 1.0d - (1.0d - siglvl^(1.0/nexp_scl))^(1.0/aexp_scl) ;Eq. 18 (anyperiod_siglvl)

end