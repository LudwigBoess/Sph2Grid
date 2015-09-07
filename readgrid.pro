;; Reads Sph2Grid output in HDF5 format.

; IDL HDF5 lib does not know about array 
; order in C and FORTRAN !
function switch_array_order, grid
    
    if (size(grid))[0] eq 3 then $
        switched_grid = transpose(grid, [2,1,0]) $
    else if (size(grid))[0] eq 4 then $
        switched_grid = transpose(grid, [0,3,2,1]) $
    else $
        switched_grid = grid

    return, switched_grid
end


function read_hdf5_grid, group_id, block

    dset_id = H5D_OPEN(group_id, block)
    dtype_id = H5D_GET_TYPE(dset_id)

    grid = H5D_READ(dset_id,dtype_id)

    H5D_CLOSE, dset_id

    return, grid
end

function readgrid, fname, block,fft=fft, Pk=Pk, head=head, $
    debug=debug, scalar=scalar
	
	if not keyword_set(fname) then begin
		print, 'grid = readgrid(fname, block,head=head, fft=fft,'
        print, '             Pk=Pk, debug=debug, scalar=scalar)'
		print, 'Returns Sph2Grid grid data '
		print, 'fname:  input filename                  '
        print, 'block:  blockname'
        print, 'fft:    return complex FFT grid instead'
        print, '        Block KVEC gives k-vector grid'
        print, 'Pk:     return binned powerspectrum'
        print, '        Block KPK gives k-vector grid'
		print, 'head:   Header structure, optional'
        print, 'scal:   return length of a vector grid'
		return, -1
	end

    if block eq 'KSCA' then begin
        block = "KVEC"
        scalar = 1
    end

	file_id = H5F_OPEN(fname)	

	; header
    head = read_hdf5_grid(file_id, "HEAD")

    if block eq "HEAD" then $
        return, head

    if arg_present(debug) then begin
        print, "File Contents :"
        nmember = H5G_GET_NUM_OBJS(file_id)

        for i=0, nmember-1 do begin
            name = H5G_GET_OBJ_NAME_BY_IDX(file_id, i)

            print, name+' <==> '+block
        end

    end

    group_id = H5G_OPEN(file_id, block) 

    ; grid
    if keyword_set(fft) then begin
        rdata = read_hdf5_grid(file_id, block+"/FFTGrid_real")
        idata = read_hdf5_grid(file_id, block+"/FFTGrid_imag")

        grid = complex(rdata, idata)

        grid = switch_array_order(grid)
    end else if keyword_set(pk) then begin
        grid = read_hdf5_grid(file_id, block+"/Pk")
    end else begin
        grid = read_hdf5_grid(group_id, "Grid")
        grid = switch_array_order(grid)
    end

    H5G_CLOSE, group_id
	H5F_CLOSE, file_id

    if keyword_set(scalar) then $
        grid = reform(sqrt(grid[0,*,*,*]^2 $
                         + grid[1,*,*,*]^2 $
                         + grid[2,*,*,*]^2))
    
	return, grid
END

