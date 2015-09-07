; plot velocity powerspectrum on a grid made by Sph2Grid
; and compare to IDL FFT. Note that differences at small k are due
; to the different data layout. IDL includes the redundant Hermitian
; part of the FFT data, which gives differences in the binning.

pro powerspectrum, fname=fname, block=block
    
    if not keyword_set(fname) then $
		fname = './grid_000'
	
	if not keyword_set(block) then $
		block = 'VEL'

    velgrid = readgrid(fname, block, head=head)

    gridsize = head.gridsize
    ngrid = head.ngrid
    cellsize = gridsize / ngrid

      ; check for empty cells
    bad = where(finite(velgrid) eq 0 ,nEmpty)
    print, 'Found empty cells :'+strn(nEmpty)

    kmin = 2*!pi/(gridsize)      ; box mode
    kmax = !pi*ngrid/gridsize    ; Nyquist mode

    kGrid = readgrid(fname, "KVEC",/scal)

    kData = readgrid(fname, block, /fft)
    kData = abs(kData[*,*,*,*])^2	; P(k) = |f^(k)|^2 for f^(k) complex
 	
	; 3D -> 1D
    lkData = sqrt(kData[0,*,*,*]^2 + kData[1,*,*,*]^2 +kData[2,*,*,*]^2)
 
    ; bin  
    good = where(kgrid le kmax)
    Pk = bin_dist(lkData[good], pos=kgrid[good], bin_pos=bin_pos, $
		nbins=20,/log,/silent)
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	Pk = readgrid(fname, "VEL", /pk)
    bin_pos = readgrid(fname, "KPK")

    ; plot
    plot, dist(2), /nodata, /xlog, /ylog, xrange=[kmin, kmax] $
		,xtitle='k = 1/l [kpc!U-1!N]' , ytitle='P(k) ' $
       , xstyle=1,ystyle=1, yrange=minmax(Pk[where(Pk ne 0)])

    oplot, bin_pos, Pk,col=color(0)
    
    oplot, bin_pos, bin_pos^(-11./3.)/1.2e15
stop
end


pro powerspectrum_compare, fname=fname
    
    if not keyword_set(fname) then fname = './grid_000'

    velgrid = readgrid(fname, 'VEL', head=head)

    gridsize = head.gridsize
    ngrid = head.ngrid
    cellsize = gridsize / ngrid

    ; check for empty cells
    bad = where(finite(velgrid) eq 0 ,nEmpty)
    print, 'Found empty cells :'+strn(nEmpty)

    kmin = 2*!pi/(gridsize)      ; box mode
    kmax = !pi*ngrid/gridsize    ; Nyquist mode

    IDLkGrid = make_kGrid(ngrid, kmin)
    kGrid = readgrid(fname, "KVEC",/scal)

    ; spectrum
    IDLkData = make_array(3,ngrid,ngrid,ngrid)
	IDLkData[0,*,*,*] = abs(FFT(velgrid[0,*,*,*],/double))^2.
    IDLkData[1,*,*,*] = abs(FFT(velgrid[1,*,*,*],/double))^2.
	IDLkData[2,*,*,*] = abs(FFT(velgrid[2,*,*,*],/double))^2.
	
    kData = readgrid(fname, "VEL", /fft)
    kData = abs(kData[*,*,*,*])^2
 
    ; bin
    lkData = sqrt(kData[0,*,*,*]^2 + kData[1,*,*,*]^2 +kData[2,*,*,*]^2)
    IDLlkData = sqrt(IDLkData[0,*,*,*]^2 + IDLkData[1,*,*,*]^2 + IDLkData[2,*,*,*]^2)

    good = where(kgrid le kmax)
    Pk = bin_dist(lkData[good], pos=kgrid[good], bin_pos=bin_pos, nbins=100,/log,/silent)
    good = where(IDLkgrid le kmax)
    idlPk = bin_dist(IDLlkData[good], pos=IDLkgrid[good], bin_pos=IDLbin_pos, nbins=100,/log,/silent)

    ; sph2grid spectrum
    sPk = readgrid(fname, "VEL", /pk)
    ksPk = readgrid(fname, "KPK")

    ; plot
    plot, dist(2), /nodata, $
       /xlog, /ylog, xrange=[kmin, kmax], yrange=[1d-15, 1d-5] $
       ,xtitle='k = 1/l [kpc!U-1!N]' $
       , ytitle='P(k) ' $
       , xstyle=1,ystyle=1

    oplot, bin_pos, Pk,col=color(0), psym=4, symsize=0.3
    oplot, IDLbin_pos, IDLPk,col=color(1), psym=5, symsize=0.3
    oplot, ksPk, sPk,col=color(2), psym=6, symsize=0.3
    
    oplot, bin_pos, bin_pos^(-11./3.)/1e16


    stop
	
	return
end

;Generate k values of an 3D grid for FFTs
function make_kGrid, nmesh, kmin
	kmag = make_array(nmesh,nmesh,nmesh,/double)

	for i=0,nmesh-1 do $
	for j=0,nmesh-1 do $
	for k=0,nmesh-1 do begin
;Define conjugated indizes of the grid
			if i ne 0 then iconj = nmesh - i $
					  else iconj = 0
			if j ne 0 then jconj = nmesh - j $
					  else jconj = 0
			if k ne 0 then kconj = nmesh - k $
					  else kconj = 0
;Define grid
			if i LE nmesh/2. then kx = i * kmin $
							 else kx = (iconj) * kmin
	
			if j LE nmesh/2. then ky = j * kmin $
							 else ky = (jconj) * kmin
	
			if k LE nmesh/2. then kz = k * kmin $
							 else kz = (kconj) * kmin

		kmag[i,j,k]  = length([kx,ky,kz])
        
	end

	return,kmag
end

FUNCTION bin_dist, array, pos=a, nbins=nbins, binsize=binsize, bin_pos=bin_pos ,log=log, counts=counts,maxmin=maxmin, std_dev=std_dev,median=median,silent=silent,rms=rms,total=total,var=var,sqsum=sqsum

	IF NOT keyword_set(array) THEN BEGIN
		print, 'bin_dist(array, pos=r, nbins=nbins, binsize=binsize, bin_pos=bin_pos,counts=counts)'
		print, '---------------------------------------------------------------------'
		print, 'Bins a 1D array using averages.'
		print, 'pos	   -  Distances, if array is not sorted'
		print, 'log	   -	Bin logarithmically'
		print, 'nbins	-	Number of bins to use'
		print, 'binsize- 	Size of the bins'
		print, 'bin_pos-	Array with start of bin positions'
		print, 'counts	-  Array with N_elements for every bin'
		print, 'maxmin	-  [2,nbins] array with max and min in each bin'
		print, 'std_dev-	std deviation in each bin'
		print, '---------------------------------------------------------------------'
		print, 'If only one of nbins/binsize is given the other one is calculated'
		print, 'assuming the full range of the array'
		print, 'IF both are given, the array is binned outwards up to nbins*binsize' 
		return, -1
	END

	if keyword_set(a) THEN $
	    pos = a

	IF NOT keyword_set(silent)THEN silent = 0 ELSE silent = 1
	IF NOT keyword_set(pos) 	THEN pos	= findgen(n_elements(array))  
	
	offset	= min(pos)
	;print, 'Offset is : '+strn(offset)

	IF keyword_set(log) THEN BEGIN
		IF NOT SILENT THEN print, 'Logarithmic Binning'
		IF min(pos) EQ 0. THEN BEGIN 
			IF NOT SILENT THEN print, 'Neglecting Datapoint at r=0 !'
			pos	= pos[where(pos NE 0)]
			offset	= min(pos)			
			IF NOT SILENT THEN print, 'NEW Offset is : '+strn(offset)
		END
		IF keyword_set(nbins) AND keyword_set(binsize) THEN $
            IF NOT SILENT THEN $
                print, 'Maximum radius : '+strn(10.^(nbins*binsize))
		IF keyword_set(nbins) AND NOT keyword_set(binsize) THEN $
            binsize	= alog10(max(pos)/min(pos))/nbins
		IF NOT keyword_set(nbins) AND keyword_set(binsize) THEN $
            nbins	= ceil(alog10(max(pos)/min(pos))/nbins)

		bin_pos	= offset*10^(binsize*lindgen(nbins))
			;10.^(findgen(nbins)*binsize)+offset-1.

	END ELSE BEGIN
		IF 	 keyword_set(nbins)	AND 	 	keyword_set(binsize) THEN IF NOT SILENT THEN print, 'Maximum radius : '+strn(nbins*binsize)
		IF 	 keyword_set(nbins)	AND NOT 	keyword_set(binsize) THEN binsize	= (max(pos)-min(pos))/nbins
		IF NOT keyword_set(nbins) AND 		keyword_set(binsize) THEN nbins		= ceil(max(pos)/binsize)
		bin_pos	= findgen(nbins)*binsize +offset + 0.5*binsize
	END
	
	profile	= make_array(nbins, value=0.,/double)
	median	= make_array(nbins, value=0.,/double)
	counts	= make_array(nbins, value=0.,/double)
	maxmin	= make_array(2,nbins,value=0.)
	std_dev	= make_array(nbins,value=0.)
	rms 	= make_array(nbins,value=0.)
	total	= make_array(nbins,value=0.,/double)
	var		= make_array(nbins,value=0.,/double)
	sqsum	= make_array(nbins,value=0.,/double)
	

	FOR i=0L, nbins-1 DO BEGIN
		IF keyword_set(log) THEN BEGIN 
			good	= where(pos GT 10^( i	*binsize)*offset AND   		$
							  pos LT 10^((i+1)*binsize)*offset,count) 	
		END ELSE BEGIN 
			good	= where(pos GT  i	  		*binsize	+offset  AND 			$
			              pos LT (i+1)		*binsize	+offset,count)
		END
		IF count NE 0 THEN BEGIN
			profile[i]	= 	mean(array[good],/double)
			total[i]	= TOTAL(array[good],/double)
			maxmin[*,i]	= [max(array[good]),min(array[good])]
			median[i]	= MEDIAN(array[good],/double)
			IF count NE 1 THEN BEGIN
				std_dev[i] 	=  STDDEV(array[good],/double)
				rms[i]		= sqrt(mean(array[good]^2,/double))
				var[i]		= variance(array[good],/double)
				sqsum[i]			= TOTAL(array[good]^2,/double)
			END
		END
		counts[i]	= count
	END

return, profile
END
