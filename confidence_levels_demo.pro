PRO MODELJ, x, a, f, pder  ;This function implements the noise model of Equation 26

  ax = a[0]*x^a[1]  ;Exponential component

  IF (a[0] LT 0) OR (a[2] LT 0) OR (a[3] LT 0) OR (a[4] LT 0) OR (a[5] LT 0) THEN $
  	f = 0*x - 1d10 $ ;discard negative amplitude solutions
  ELSE BEGIN
  	IF (a[3] EQ 0.0) or (a[4] EQ 0.0) THEN $  ;IF cases kappa = 0 or rho = 0 THEN Kappa function = 0
  		bx = 0.0 $
  	ELSE BEGIN
  		IF a[3] GT 1d2 THEN a[3] = 1d2	;Kappa function does not vary significantly if Kappa > 100 (tends towards a Gaussian). This helps convergence.
  		bx = a[2]*(1.0d + (x^2.0d)/(a[3]*a[4]^2.0d))^(-(a[3]+1.0d)/2.0d)  ;Kappa function (kappa = a[3])
  	ENDELSE
  	f = ax + bx + a[5]  ;Exponential + Kappa function + constant
  ENDELSE

  IF N_PARAMS() GE 4 THEN $ ;Partial derivatives (Mathematica !)
    pder = [[ax/a[0]], $
			[ax*ALOG(x)], $
			[bx/a[2]], $
			[bx*((a[3]+1.0d)*x^2.0d - a[3]*(x^2.0d + a[3]*a[4]^2.0d)*ALOG(1.0d + (x^2.0d)/(a[3]*a[4]^2.0d)))/(2.0d*a[3]*(x^2.0d + a[3]*a[4]^2.0d))], $
			[(a[2]*a[3]*a[4]*(1.0d + a[3])*(x^2.0d)*(1.0d + (x^2.0d)/(a[3]*a[4]^2.0d))^((1.0d - a[3])/2.0d))/((x^2.0d + a[3]*a[4]^2.0d)^2.0d)], $
			[REPLICATE(1.0d, N_ELEMENTS(x))]]

END

infile = 'cas_1_light_curve_335_demo.save'
outfile = 'figure_wavelet_' + STRMID(infile, 0, STRPOS(infile, '.save')) + '.eps'

restore, infile	;provided save file contains two arrays: light_curve (the data) and time (in hours)

dt = time[1] - time[0]

n = N_ELEMENTS(light_curve)

avg_light_curve = TOTAL(light_curve)/n
light_curve = light_curve - avg_light_curve ;Subtract the mean before any wavelet or Fourier analysis

vari = (MOMENT(light_curve))[1]
light_curve = light_curve / SQRT(vari) ;Normalize the light curve to its standard deviation


;-------------Wavelet spectrum ---------------------------------------

siglvl = 0.95 ;95% global significance level, i.e. taking into account the total number of degrees of freedom in the Fourier or wavelet spectra

;Empirical coefficients for the Morlet wavelet, see the Torrence & Compo paper. Will be different for other wavelets
mother = 'Morlet'
k0 = 6.0d
morlet_fourier_factor = (4.0*!DPI)/(k0 + SQRT(2.0+k0^2.0d))
morlet_coi = morlet_fourier_factor/SQRT(2)
gamma_m = 2.32d ;empirical factor for time averaging Morlet
morlet_dj0 = 0.6d ;empirical factor for scale averaging Morlet

s0 = 2.0*dt
noct = ALOG(n)/ALOG(2.0) - 1
dj = 1/8.0
j1 = FIX(noct/dj)

wave = WAVELET(light_curve, dt, PERIOD=period, SCALE=scale, S0=s0, $
                PAD=1, COI=coi, DJ=dj, J=j1, MOTHER=mother)

j = N_ELEMENTS(scale) - 1

power = (ABS(wave))^2.0d

global_ws = TOTAL(power, 1, /DOUBLE, /NAN)/n 	;time-averaged wavelet spectrum (Global Wavelet Spectrum in TC98)

;-----------------------------------Fourier Spectrum & Fit with EQ 26 --------------------------

model_function = 'MODELJ'
params = [0.05, -2.0, 20, 10.0, 0.1, 1d-4]	 ;Initial guess

;Used to select to various components of the fit
kappa_mask = [0, 0, 1, 1, 1, 0]
powerlaw_mask = [1, 1, 0, 0, 0, 0]
background_mask = [1, 1, 1, 1, 1, 1]
noise_mask = [0, 0, 0, 0, 0, 1]

apodized = light_curve * HANNING(n, ALPHA=0.5)  ;Hann apodized version of the light_curve

ff = FFT(apodized, /DOUBLE)

;Normalize the Fourier power spectrum with the convention used by Torrence & Compo so that the
;Fourier and time averaged-wavelet spectra line-up
vari = (MOMENT(apodized))[1]
psdt = n*DOUBLE(ff*CONJ(ff))/vari

psdt = psdt[1:n/2]
nu_psdt = DINDGEN(n)/(dt*n)
nu_psdt = nu_psdt[1:n/2]	;Frequencies at which the FFT is computed

;The power spectrum is fitted using curvefit. In order for the least squares fit to converge, the
;residuals (differences between data points and fit) must be comparable at each frequency,
;This is achieved by normalizing the squared residuals by the variance of the data.
;At each frequency nu, the Fourier spectrum is distributed as 0.5*sigma(nu)*chi2_2 (degree 2 chi2, Eq. 11)
;Its standard deviation is thus (0.5*sigma(nu))*2=sigma(nu). The standard deviation of the power
;spectrum is thus equal to its mean at each frequency. Therefore we need to weight the fit by the result
;of the fit squared (the mean power squared as a function of frequency) ... But that is OK because only
;a first order estimate of the mean power is needed. We can simply use the time-averaged wavelet spectrum as a proxy.
weights = 1.0/(INTERPOL(global_ws, period, 1/nu_psdt)^2.0d)
fit_res = CURVEFIT(nu_psdt, psdt, weights, params, FUNCTION_NAME=model_function, /DOUBLE, ITMAX=100, TOL=1d-4)

CALL_PROCEDURE, model_function, nu_psdt, params*powerlaw_mask, powerlaw_fit  ;these lines generate the various components of the resulting fit
CALL_PROCEDURE, model_function, nu_psdt, params*kappa_mask, kappa_fit
CALL_PROCEDURE, model_function, nu_psdt, params*noise_mask, noise_fit
CALL_PROCEDURE, model_function, nu_psdt, params*background_mask, background_fit

CALL_PROCEDURE, model_function, 1.0/period, params*background_mask, background_fit_period

;-------------Confidence Levels ---------------------------------------

jcoi = ALOG(coi/morlet_fourier_factor/s0)/ALOG(2)/dj ;scale index (j) corresponding to each period of the COI

;-- Global Fourier confidence level (Eq. 12)
anyperiod_fourier_tresh = -ALOG(1.0d - siglvl^(1.0/(n/2.0d)))

;-- Local confidence level for the wavelet spectrum (i.e. as in TC98, but with a custom background model (cf. GWS keyword))
powerlaw_signif = WAVE_SIGNIF(light_curve, dt, scale,0, $
			  GWS=background_fit_period, $
              SIGLVL=siglvl,MOTHER=mother)
powerlaw_signif = REBIN(TRANSPOSE(powerlaw_signif), n, J+1)  ; expand powerlaw_signif to the size of the wavelet spectrum

powerlaw_signif = power/powerlaw_signif   ; where powerlaw_signif > 1, there is power above the local confidence level

 ;-- Global confidence level for the wavelet spectrum

Nout = TOTAL(jcoi>0)  ;Number of points outside the COI in the wavelet spectrum
anydtperiod_siglvl = GLOBAL_TO_LOCAL_SIGLVL(siglvl, Mother, Nout, dj) ;Local confidence level in the wavelet spectrum for the desired global confidence level

anydtperiod_powerlaw_signif = WAVE_SIGNIF(light_curve,dt,scale,0.0, $
			 gws=background_fit_period, $
             SIGLVL=anydtperiod_siglvl,MOTHER=mother)
anydtperiod_powerlaw_signif = REBIN(TRANSPOSE(anydtperiod_powerlaw_signif),n,J+1)  ; expand to the size of the wavelet spectrum

anydtperiod_powerlaw_signif = power/anydtperiod_powerlaw_signif   ; where anydtperiod_powerlaw_signif > 1, there is power above the global confidence level

;-- Local confidence level for the time-averaged wavelet spectrum (i.e. as in TC98, but with a custom background model (cf. GWS keyword))
dof = (n - scale/dt)>0
gws_powerlaw_signif = WAVE_SIGNIF(light_curve, dt, scale, 1, $
		   SIGLVL=siglvl, GWS=background_fit_period, $
           DOF=dof, MOTHER=mother)

;-- Global confidence level for the time-averaged wavelet spectrum

Sout = MAX(jcoi)	;Number of points outside the COI in the time-averaged wavelet spectrum
anyperiod_siglvl = GLOBAL_TO_LOCAL_SIGLVL(siglvl, Mother, Sout, dj, /GWS) ;Local confidence level in the time-averaged wavelet spectrum (/GWS) for the desired global confidence level

dof = (n - scale/dt)>0
anyperiod_gws_powerlaw_signif = WAVE_SIGNIF(light_curve, dt, scale, 1, $
			SIGLVL=anyperiod_siglvl, GWS=background_fit_period, $
			DOF=dof, MOTHER=mother)


whitened_power = power / REBIN(TRANSPOSE(background_fit_period), n, J+1) ;For each time step, normalize (whiten) the wavelet power by the noise model


;--------------------------- Plotting --------------------------------------------------

dt_unit = 'hours'
freq_unit = '!9' + STRING(109b) + '!3Hz'

power_label = 'Power (!9' + STRING(115b) + '!S!3' + STRING(178b) + '!R!D0!N)'

ps_xsize = 20.0
ps_ysize = 9.7

ps_yscale = 16.0/ps_ysize

wvlplthgt = 0.34*ps_yscale

tick_len = 0.05

thin_line = 1
medium_line = 2
thick_line = 3

cmax = 240
LOAD_WAVELET_COLORS

red_mean_color = !blu
red_local_wavelet_color = !orange
red_global_wavelet_color = !yellow
red_global_fourier_color = !yellow

powerlaw_mean_color = !red
powerlaw_local_wavelet_color = !orange
powerlaw_global_wavelet_color = !yellow
powerlaw_global_fourier_color = !light_grey

fourier_mean_color = !grey

mydevice = !d.name

SET_PLOT, 'ps'
DEVICE, FILE = outfile, /ENCAPSULATED, /COLOR, $
		BITS=8, /HELVETICA, XS=ps_xsize, YS=ps_ysize, /TT_FONT

	;--- Plot time series
pos1 = [0.05,1-0.15*ps_yscale,0.5,1.0]
xrange_time = [0, n*dt]
yrange = [-1.0, 1.0]*CEIL(MAX(ABS(light_curve)))

PLOT,time,light_curve, XRANGE=xrange_time, YRANGE=yrange, POSITION=pos1, $
	YTITLE='Intensity (!9' + STRING(115b) + '!3!D0!N)', $
	XS=1+4, YS=1+8, THICK=thin_line, COLOR=!black, /NODATA, YMINOR=1, /NOERASE, $
	YTICKLEN=-tick_len/ps_xsize/(pos1[2] - pos1[0]), FONT=1

OPLOT, time, light_curve, COLOR=!grey, THICK=medium_line

minval = 0.0
maxval = 0.9

minim = MIN(whitened_power, MAX = maxim, /NAN)
maxim = maxim<100
histo = HISTOGRAM(whitened_power, MIN = minim, MAX = maxim, NBINS=100, LOCATIONS = x)
tot = TOTAL(ABS(histo*x), /DOUBLE, /CUMULATIVE)
tot = tot/tot[N_ELEMENTS(tot)-1]

good = WHERE((tot GE minval) AND (tot LE maxval), count)
IF count GT 0 THEN $
	mini = MIN(x[good], MAX = maxi) $
ELSE mini = MIN(whitened_power, MAX = maxi)

IF mini EQ maxi THEN mini = MIN(whitened_power, MAX = maxi)

max_power = maxi
noct = 10.0^FLOOR(ALOG10(maxi)-0.0)
max_power = noct*FLOOR(max_power/noct)

;--- Contour plot wavelet power spectrum
yrange = [1.0/MIN(period)/3600.0, 1.0/MAX(period)/3600.0]*1d6   ; Hz
nlevels = 256
disp_power = whitened_power
mini = 0
maxi = max_power

levels = mini + (maxi - mini)*DINDGEN(nlevels)/nlevels
colors = cmax*DINDGEN(nlevels)/nlevels
period2 = FLOOR(ALOG(1d6/period/3600.0)/ALOG(10))   ; integer powers of 2 in period
log_ytickv = period2[UNIQ(period2)]
ytickv = 10.0^log_ytickv  	; unique powers of 2
gd = WHERE((ytickv GE MIN(yrange)) AND (ytickv LE MAX(yrange)))
ytickv = ytickv[gd]
log_ytickv = log_ytickv[gd]
nyticks = N_ELEMENTS(ytickv)
power_tickname = '10!U'+ STRTRIM(STRING(log_ytickv, format='(I2)'), 2) + '!N'

xtickv = 20*INDGEN(10)
gd = WHERE((xtickv GE MIN(xrange_time)) AND (xtickv LE MAX(xrange_time)))
xtickv = xtickv[gd]

xtitle = 'Time (' + dt_unit + ')'

pos2 = [pos1(0),0.07*ps_yscale,pos1(2),0.07*ps_yscale + wvlplthgt]

IF (J+1) LT 512 THEN wvl_ydisp = CEIL(512.0/(J+1))*(J+1) ELSE wvl_ydisp = J+1
IF n LT 512 THEN wvl_xdisp = CEIL(512.0/n)*n ELSE wvl_xdisp = n

disp_power = REBIN(disp_power, wvl_xdisp, wvl_ydisp)

IF n GT 1024 THEN wvl_xdisp = 512.0

disp_power = CONGRID(disp_power, wvl_xdisp, wvl_ydisp)

disp_mask = REPLICATE(1.0, wvl_xdisp, wvl_ydisp)
rscale = REBIN(scale, wvl_ydisp)
FOR i=0, wvl_ydisp-1 DO BEGIN
	disp_mask[0:(rscale[i]/(morlet_coi*n/wvl_xdisp)/dt)<(wvl_xdisp/2), i] = 4.0
	disp_mask[wvl_xdisp-((rscale[i]/(morlet_coi*n/wvl_xdisp)/dt)<(wvl_xdisp/2)):*, i] = 4.0
ENDFOR

disp_power = BYTE(ROUND(cmax*((disp_mask*disp_power)<maxi - mini)/(maxi - mini)))

PLOT,[1, 2], [1, 2],/NOERASE,POSITION=pos2, $
         XRANGE=xrange_time,YRANGE=yrange,/YTYPE, $
         XTICKS=1, yticks=1, XS=1+4, YS=1+4, COLOR=!black, /NODATA

TV, disp_power, !X.WINDOW[0], !Y.WINDOW[0], $
		    XSIZE = !X.WINDOW[1] - !X.WINDOW[0], $
		    YSIZE = !Y.WINDOW[1] - !Y.WINDOW[0], /NORMAL

PLOT, [1, 2], [1, 2],/NOERASE,POSITION=pos2, $
         XRANGE=xrange_time,YRANGE=yrange,/YTYPE, $
         YTICKS=N_ELEMENTS(ytickv)-1,YTICKV=ytickv, $
         YTICKNAME=power_tickname, $
         YMINOR=9, $
         XMINOR=2, $
         XTICKLEN=-tick_len/ps_ysize/(pos2[3] - pos2[1]), $
         YTICKLEN=-tick_len/ps_xsize/(pos2[2] - pos2[0]), $
         XTITLE=xtitle,YTITLE=' ', XTICKFORMAT = '(I3)', $
         XTICKS=N_ELEMENTS(xtickv)-1, $
         XTICKV=xtickv, /XS, /YS, COLOR=!black, /NODATA, FONT=1

axis, yaxis=0, pos2[0], pos2[1], /YLOG, COLOR=!black, $
	YTICKS=1, YTICKNAME=['   ', '   '], YTITLE='Frequency (' + freq_unit + ')', FONT=1, /NORMAL, $
	YTICKLAYOUT=1

CONTOUR,powerlaw_signif,time, 1d6/period/3600.0,/OVERPLOT,LEVEL=1,THICK=medium_line, $
        C_COLORS=powerlaw_local_wavelet_color

CONTOUR,anydtperiod_powerlaw_signif,time,1d6/period/3600.0,/OVERPLOT,LEVEL=1,THICK=medium_line, $
		C_COLORS=powerlaw_global_wavelet_color

;Produces the two right-side plots
FOR ilog = 0, 1 DO BEGIN

	IF ilog EQ 0 THEN xtitle = 'Normalized power (!9' + STRING(115b) + '!3)' ELSE xtitle = power_label

	IF ilog EQ 0 THEN BEGIN
		norm = background_fit
		norm_period = background_fit_period
	ENDIF ELSE BEGIN
		norm = 1.0
		norm_period = 1.0
	ENDELSE

	IF ilog EQ 1 THEN xmin = 1e-4 ELSE xmin = 1e-2
	xmax = 10.0^CEIL(ALOG10(MAX([anyperiod_fourier_tresh*background_fit/norm, psdt/norm])))
	IF ilog EQ 0 THEN xmax = xmax*10.0
	IF ilog EQ 1 THEN xmax = 1e4 ELSE xmax = 1e3

	xrange = [xmin, xmax]

	XTICKS = (ALOG10(xrange[1]) - ALOG10(xrange[0]))>1
	xtickname = STRARR(XTICKS+1)
	xtickname = '10!U'+ STRTRIM(STRING(ALOG10(xrange[0]) + INDGEN(XTICKS+2), format='(I2)'), 2) + '!N'

	;--- Plot Fourier & time-averaged wavelet spectrum
	pos3 = [pos2[2]+ilog*0.25+0.04,pos2[1],pos2[2]+ilog*0.25+0.04+0.21,pos2[3]]

	PLOT, psdt,(1/nu_psdt)/3600.0,/NOERASE, POSITION=pos3, COLOR=!black, $
			/XLOG, $
			/YLOG, $
		 	XSTYLE=8+1,YSTYLE=8+1, $
			XRANGE=xrange, $
		    YRANGE=yrange, /YTYPE, $
		    XTICKS=XTICKS,XMINOR=1, PSYM=10, $
		    XTICKNAME=xtickname, $
		    YTICKS=N_ELEMENTS(ytickv)-1,$
		    YTICKV=ytickv,$
		    YMINOR=9, $
		    YTICKNAME=power_tickname, $
			XTICKLEN=-tick_len/ps_ysize/(pos3[3] - pos3[1]), $
			YTICKLEN=-tick_len/ps_xsize/(pos3[2] - pos3[0]), $
		    XTITLE=xtitle, /NODATA, /NORMAL, FONT=1

	;Histogram style Fourier spectrum
	npsdt = N_ELEMENTS(psdt)
	histo_psdt = CONGRID(psdt/norm, npsdt*2)
	period_psdt = 1/nu_psdt
	histo_period_psdt = CONGRID(period_psdt, npsdt*2, /CENTER)
	PLOTS, histo_psdt,1d6/histo_period_psdt/3600.0,COLOR=!grey,/DATA, NOCLIP=0

	;Global wavelet spectrum
	OPLOT, global_ws/norm_period, 1d6/period/3600.0,THICK=thin_line, COLOR=!black

	;Significance levels
	OPLOT, background_fit_period/norm_period, 1d6/period/3600.0, COLOR=powerlaw_mean_color, THICK=thin_line
	OPLOT, gws_powerlaw_signif/norm_period,1d6/period/3600.0, COLOR=powerlaw_local_wavelet_color, THICK=thin_line
	OPLOT, anyperiod_gws_powerlaw_signif/norm_period,1d6/period/3600.0, COLOR=powerlaw_global_wavelet_color, THICK=thin_line
	OPLOT, anyperiod_fourier_tresh*background_fit_period/norm_period, 1d6/period/3600.0, COLOR=powerlaw_global_fourier_color, THICK=medium_line

	;Background noise fit & components
	OPLOT, kappa_fit/norm,1d6/period_psdt/3600.0, COLOR=!red, THICK=thin_line, LINE=2
	OPLOT, powerlaw_fit/norm,1d6/period_psdt/3600.0, COLOR=!red, THICK=thin_line, LINE=2
	OPLOT, noise_fit/norm,1d6/period_psdt/3600.0, COLOR=!red, THICK=thin_line, LINE=2

ENDFOR

DEVICE, /CLOSE
SET_PLOT, mydevice

END