pro load_wavelet_colors

	cmax = 240.0
	red_idx = [32 + (255.0 - 32)*indgen(cmax+1)/cmax, replicate(0.0, 255 - cmax)]
	gre_idx = red_idx
	blu_idx = red_idx

	defsysv, '!red', 255
	defsysv, '!orange', 254
	defsysv, '!yellow', 253
	defsysv, '!dark_blu', 252
	defsysv, '!blu', 251
	defsysv, '!bright_blu', 250
	defsysv, '!dark_gre', 249
	defsysv, '!gre', 248
	defsysv, '!bright_gre', 247
	defsysv, '!magenta', 246
	defsysv, '!cyan', 245
	defsysv, '!grey', 244
	defsysv, '!light_grey', 243
	defsysv, '!black', 242
	defsysv, '!white', 241

	red_idx[!red] = 255
	gre_idx[!red] = 0
	blu_idx[!red] = 0

	red_idx[!orange] = 255
	gre_idx[!orange] = 128
	blu_idx[!orange] = 64

	red_idx[!yellow] = 255
	gre_idx[!yellow] = 192
	blu_idx[!yellow] = 128

	red_idx[!dark_blu] = 0
	gre_idx[!dark_blu] = 0
	blu_idx[!dark_blu] = 128

	red_idx[!blu] = 64
	gre_idx[!blu] = 64
	blu_idx[!blu] = 255

	red_idx[!bright_blu] = 64
	gre_idx[!bright_blu] = 128
	blu_idx[!bright_blu] = 255

	red_idx[!dark_gre] = 55
	gre_idx[!dark_gre] = 100
	blu_idx[!dark_gre] = 0

	red_idx[!gre] = 82
	gre_idx[!gre] = 150
	blu_idx[!gre] = 0

	red_idx[!bright_gre] = 110
	gre_idx[!bright_gre] = 200
	blu_idx[!bright_gre] = 0

	red_idx[!magenta] = 255
	gre_idx[!magenta] = 64
	blu_idx[!magenta] = 255

	red_idx[!cyan] = 0
	gre_idx[!cyan] = 200
	blu_idx[!cyan] = 200

	red_idx[!grey] = 80
	gre_idx[!grey] = 80
	blu_idx[!grey] = 80

	red_idx[!light_grey] = 200
	gre_idx[!light_grey] = 200
	blu_idx[!light_grey] = 200

	red_idx[!black] = 0
	gre_idx[!black] = 0
	blu_idx[!black] = 0

	red_idx[!white] = 255
	gre_idx[!white] = 255
	blu_idx[!white] = 255

	tvlct, red_idx, gre_idx, blu_idx

end