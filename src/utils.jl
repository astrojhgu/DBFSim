module utils
    using WCS
    using AstroLib
    
    
    export meshgrid

    function meshgrid(x, y)
        X = [i for i in x, j in 1:length(y)]
        Y = [j for i in 1:length(x), j in y]
        return X, Y
    end

    function azel2xyz(az, el)
        az=deg2rad(90.0-az)
        el=deg2rad(el)
        ce=cos(el)
        se=sin(el)
        ca=cos(az)
        sa=sin(az)
        x=ce*ca
        y=ce*sa
        z=se
        [x,y,z]
    end
    


    function sph_grid(naz::Int, nel::Int)
        az_grid, el_grid=meshgrid(range(0.0, 360.0, length=naz), range(0.0, 90.0, length=nel))
        #azs=reshape(az_grid, naz*nel)
        #els=reshape(el_grid, naz*nel)

        #xyz=hcat(azel2xyz.(azs, els)...)'
        #azs,els 

        xyz=map(zip(az_grid, el_grid)) do(az, el)
            azel2xyz(az, el)
        end

        x=map(xyz) do xyz
            xyz[1]
        end
        y=map(xyz) do xyz
            xyz[2]
        end
        z=map(xyz) do xyz
            xyz[3]
        end

        az_grid, el_grid, x, y, z
    end

    function small_field_grid(ha, dec, lat, field_size_deg, npix)
        dx=field_size_deg/npix
        wcs = WCSTransform(2;
            cdelt = [-dx, dx],
            ctype = ["RA---ARC", "DEC--ARC"],
            crpix = [0., 0.],
            crval = [ha, dec])

        az=zeros(npix, npix)
        el=zeros(npix, npix)

        for i in 1:npix
            for j in 1:npix
                x=j-npix/2
                y=i-npix/2
                sky_point=pix_to_world(wcs, [x, y])
                el1, az1=hadec2altaz(sky_point... , lat)
                el[i,j]=el1
                az[i,j]=az1
            end
        end
        #xyz=hcat(azel2xyz.(reshape(az, npix^2), reshape(el, npix^2))...)'
        az, el
    end


    function gen_faces(az, el, g)
        naz=size(az, 1)
        nel=size(az, 2)
        indices=Vector{Int}()
        for i in 1:naz-1
            for j in 1:nel-1
                k11=(j-1)*naz+i
                k12=j*naz+i
                k21=(j-1)*naz+i+1
                k22=j*naz+i+1
                push!(indices, k11, k12, k21, k12, k21, k22)
            end
        end

        i=naz
        for j in 1:nel-1
            k11=(j-1)*naz+naz
            k12=j*naz+naz
            k21=(j-1)*naz+1
            k22=j*naz+1
            push!(indices, k11, k12, k21, k12, k21, k22)
        end

        xyz=map(zip(az, el, g)) do(az, el, g)
            azel2xyz(az, el)*g
        end

        xyz=reshape(xyz, naz*nel)

        x=map(xyz) do xyz
            xyz[1]
        end
        y=map(xyz) do xyz
            xyz[2]
        end
        z=map(xyz) do xyz
            xyz[3]
        end


        c=reshape(g, naz*nel)
        c./=maximum(c)
        x,y,z,indices,c

    end

    
end