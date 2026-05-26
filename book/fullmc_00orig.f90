program test_fullmc
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    character(len=256) :: line
    integer :: nphoton
    integer :: nx, ny, nz
    real(dp) :: dx, dy, dz
    integer, allocatable :: nxs(:)
    real(dp), allocatable :: dxs(:)
    real(dp), allocatable :: xarr(:), yarr(:), zarr(:)
    real(dp), allocatable :: kext(:,:,:), ksca(:,:,:)
    real(dp), allocatable :: kabs(:,:,:), gparam(:,:,:)
    real(dp), allocatable :: galb(:,:)
    real(dp), allocatable :: bplnk(:,:,:), bgrnd(:,:)
    integer :: source, swlw, transfer_mode
    real(dp) :: solmu, solphi
    real(dp) :: viewmu, viewphi
    integer :: scaord
    real(dp) :: init_weight
    real(dp), allocatable :: dirsol(:), dirview(:)
    integer :: iphoton, it, it1, itmax, it1max
    real(dp) :: maxx, maxy, maxz
    real(dp), allocatable :: rloc(:), rdir(:)
    integer, allocatable :: rind(:), rind2(:), rindsrc(:)
    real(dp), allocatable :: r1loc(:), r1dir(:)
    integer, allocatable :: r1ind(:), r1ind2(:)
    real(dp) :: weight, ptau, rnd, phradval, phrad1val
    real(dp) :: ksca_grid, kabs_grid, kext_grid, g_grid
    integer :: ix, iy, iz, ia
    integer :: icase, isign, idi
    integer, allocatable :: rdir_sign(:)
    integer, allocatable :: r1dir_sign(:)
    real(dp), allocatable :: xbnd(:), ybnd(:), zbnd(:), dc(:), rdloc(:)
    real(dp) :: notindsrc
    real(dp) :: rdist
    real(dp) :: weight_absorbed, new_weight
    real(dp) :: weight_min, weight_rr
    real(dp) :: g, s, mu, phi, sint, ux, uy, uz, denom
    real(dp), allocatable :: outradirr(:,:,:,:,:,:), outradconv(:,:,:,:), outradimg(:,:,:)
    real(dp) :: pi = acos(-1.0_dp)
    integer :: i
    integer :: iuconf, iutraj
    integer :: seed_size
    integer, allocatable :: seed(:)
    integer :: debug

    if (command_argument_count() >= 1) then
        call get_command_argument(1, line)
    else
        read(*,'(A)') line
    end if

    iuconf = 41
    ! open(unit=iuconf, file='out/test_fullmc/config.txt', status='old', action='read')
    open(unit=iuconf, file=trim(line)//'/config.txt', status='old', action='read')
    read(iuconf,*) nx, ny, nz
    read(iuconf,*) dx, dy, dz
    read(iuconf,*) source, swlw
    read(iuconf,*) transfer_mode
    read(iuconf,*) solmu, solphi
    read(iuconf,*) viewmu, viewphi
    read(iuconf,*) nphoton
    read(iuconf,*) debug

    write(*,*) 'Configuration: nx=', nx, ' ny=', ny, ' nz=', nz
    write(*,*) '               dx=', dx, ' dy=', dy, ' dz=', dz
    write(*,*) '               nphoton=', nphoton
    write(*,*) '               source=', source
    write(*,*) '               swlw=', swlw
    write(*,*) '               transfer_mode=', transfer_mode
    ! write(*,*) '               debug=', debug

    ! allocate(xarr(0:nx-1), yarr(0:ny-1), zarr(0:nz-1))
    allocate(xarr(0:nx), yarr(0:ny), zarr(0:nz))
    ! xarr = linspace(0, nx*dx, nx) in Python -> use same max value (nx*dx)
    do i = 0, nx-1
        xarr(i) = real(i, dp) * dx
    end do
    xarr(nx) = real(nx,dp) * dx
    do i = 0, ny-1
        yarr(i) = real(i, dp) * dy
    end do
    yarr(ny) = real(ny,dp) * dy
    do i = 0, nz-1
        zarr(i) = real(i, dp) * dz
    end do
    zarr(nz) = real(nz,dp) * dz

    allocate(kext(0:nx-1,0:ny-1,0:nz-1))
    allocate(ksca(0:nx-1,0:ny-1,0:nz-1))
    allocate(kabs(0:nx-1,0:ny-1,0:nz-1))
    allocate(gparam(0:nx-1,0:ny-1,0:nz-1))
    allocate(galb(0:nx-1,0:ny-1))
    allocate(bplnk(0:nx-1,0:ny-1,0:nz-1))
    allocate(bgrnd(0:nx-1,0:ny-1))

    do ix = 0, nx-1
        do iy = 0, ny-1
            do iz = 0, nz-1
                read(iuconf,*) kext(ix,iy,iz), kabs(ix,iy,iz), gparam(ix,iy,iz), bplnk(ix,iy,iz)
                ksca(ix,iy,iz) = kext(ix,iy,iz) - kabs(ix,iy,iz)
                ! write(*,*) 'Read cell ', ix, iy, iz, ': kext=', kext(ix,iy,iz), ' kabs=', kabs(ix,iy,iz), ' g=', gparam(ix,iy,iz)
            end do
            read(iuconf,*) galb(ix,iy), bgrnd(ix,iy)
            ! write(*,*) 'Read albedo ', ix, iy, ': galb=', galb(ix,iy)
        end do
    end do

    init_weight = 1.0_dp

    ! weight_min = 0.1_dp
    ! weight_rr = 0.1_dp
    weight_min = 0.5_dp
    weight_rr = 0.5_dp

    allocate(dirsol(0:2))
    dirsol(0) = sqrt(1.0_dp - solmu**2) * sin(solphi)
    dirsol(1) = sqrt(1.0_dp - solmu**2) * cos(solphi)
    dirsol(2) = -solmu
    ! dirsol = (/ 0.5_dp, 0.0_dp, -0.5_dp /)
    dirsol = dirsol / sqrt(sum(dirsol**2))

    allocate(dirview(0:2))
    dirview(0) = sqrt(1.0_dp - viewmu**2) * sin(viewphi)
    dirview(1) = sqrt(1.0_dp - viewmu**2) * cos(viewphi)
    dirview(2) = -viewmu
    dirview = dirview / sqrt(sum(dirview**2))

    ! RNG seeding for reproducibility (seed value ~12346)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345 + (/ (i-1, i=1,seed_size) /)
    call random_seed(put=seed)

    maxx = real(nx,dp) * dx
    maxy = real(ny,dp) * dy
    maxz = real(nz,dp) * dz

    allocate(nxs(0:2))
    nxs(0) = nx
    nxs(1) = ny
    nxs(2) = nz

    allocate(dxs(0:2))
    dxs(0) = dx
    dxs(1) = dy
    dxs(2) = dz

    allocate(rloc(0:2), rdir(0:2), rind(0:2), rind2(0:2), rindsrc(0:2))
    allocate(r1loc(0:2), r1dir(0:2), r1ind(0:2), r1ind2(0:2))
    allocate(rdir_sign(0:2), r1dir_sign(0:2))
    allocate(xbnd(0:1), ybnd(0:1), zbnd(0:1), dc(0:2), rdloc(0:2))

    allocate(outradirr(0:nx-1,0:ny-1,0:nz-1,0:swlw,0:5,0:1))
    allocate(outradconv(0:nx-1,0:ny-1,0:nz-1,0:1))
    allocate(outradimg(0:nx-1,0:ny-1,0:1))
    outradirr = 0.0_dp
    outradconv = 0.0_dp
    outradimg = 0.0_dp

    itmax = 100000
    it1max = 100000

    iutraj = 42

    if (debug >= 2) then
        open (unit=iutraj, file=trim(line)//'/photon_trajectory.txt', status='replace', action='write')
        write(iutraj,*) 'ix, iy, iz, ia, iphoton, it, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), '// &
                    'rdir(0), rdir(1), rdir(2), weight, ptau'
    end if

    write(*,*) "Starting photon transport..."

    do ix = 0, nx - 1
    do iy = 0, ny - 1
    do iz = 0, merge(merge(0, nz - 1, source <= 1), 0, source <= 3)
    do ia = 0, merge(0, 5, source /= 3)
    do iphoton = 1, nphoton
        if (source <= 1) then
            ! Top of model
            ! call random_number(rnd); rloc(0) = maxx * rnd
            ! call random_number(rnd); rloc(1) = maxy * rnd
            call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
            call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
            rloc(0) = mod(rloc(0) + maxx, maxx)
            rloc(1) = mod(rloc(1) + maxy, maxy)
            rloc(2) = maxz

            rindsrc(0) = mod(int( rloc(0) / dx ), nx)
            rindsrc(1) = mod(int( rloc(1) / dy ), ny)
            rindsrc(2) = mod(int((rloc(2) - 1.0e-8_dp) / dz ), nz)
            rind(0:2) = rindsrc(0:2)
        else if (source == 2) then
            ! Volumetric source
            call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
            call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
            call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
            rloc(0) = mod(rloc(0) + maxx, maxx)
            rloc(1) = mod(rloc(1) + maxy, maxy)
            rloc(2) = mod(rloc(2) + maxz, maxz)
            
            rindsrc(0) = mod(int( rloc(0) / dx ), nx)
            rindsrc(1) = mod(int( rloc(1) / dy ), ny)
            rindsrc(2) = mod(int( rloc(2) / dz ), nz)
            rind(0:2) = rindsrc(0:2)
        else if (source == 4) then
            rloc(0) = xarr(ix) + 0.5_dp * dx
            rloc(1) = yarr(iy) + 0.5_dp * dy
            rloc(2) = maxz

            rindsrc(0) = mod(int( rloc(0) / dx ), nx)
            rindsrc(1) = mod(int( rloc(1) / dy ), ny)
            rindsrc(2) = mod(int((rloc(2) - 1.0e-8_dp) / dz ), nz)
            rind(0:2) = rindsrc(0:2)
        end if

        if (source == 0) then 
            rdir = dirsol
        else if (source == 1) then
            ! Lambertian source
            call random_number(rnd)
            mu = sqrt(rnd)
            call random_number(rnd)
            phi = 2.0_dp * pi * rnd
            sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
            rdir(0) = sint * cos(phi)
            rdir(1) = sint * sin(phi)
            rdir(2) = -mu
        else if (source == 2) then
            ! Isotropic source
            call random_number(rnd)
            mu = 2.0_dp * rnd - 1.0_dp
            call random_number(rnd)
            phi = 2.0_dp * pi * rnd
            sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
            rdir(0) = sint * cos(phi)
            rdir(1) = sint * sin(phi)
            rdir(2) = mu
        else if (source == 4) then
            rdir = dirview
        end if
        rdir = rdir / sqrt(sum(rdir**2))
        rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))
        weight = init_weight
        call random_number(rnd); ptau = -log(max(1.0e-300_dp, rnd))

        phradval = 0.0_dp
        scaord = 0

        it = 0
        do while (it < itmax)
            if (debug > 0) then
                if (debug == 1) then
                    print *, "Step ", it, " Rloc: ", rloc, " Rind: ", rind, " Rdir: ", rdir, " Weight: ", weight, " Ptau: ", ptau
                else 
                    write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'a', 0, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), &
                        rdir(0), rdir(1), rdir(2), weight, ptau
                end if
            end if
            ksca_grid = ksca(rind(0), rind(1), rind(2))
            kabs_grid = kabs(rind(0), rind(1), rind(2))
            g_grid = gparam(rind(0), rind(1), rind(2))

            xbnd(0) = xarr(rind(0))
            xbnd(1) = xarr(rind(0)+1)
            ybnd(0) = yarr(rind(1))
            ybnd(1) = yarr(rind(1)+1)
            zbnd(0) = zarr(rind(2))
            zbnd(1) = zarr(rind(2)+1)

            dc(0) = xbnd(rdir_sign(0)) - rloc(0)
            dc(1) = ybnd(rdir_sign(1)) - rloc(1)
            dc(2) = zbnd(rdir_sign(2)) - rloc(2)

            icase = 2
            icase = merge(1, icase, abs(dc(1) * rdir(icase)) < abs(dc(icase) * rdir(1)))
            icase = merge(0, icase, abs(dc(0) * rdir(icase)) < abs(dc(icase) * rdir(0)))

            if (debug == 1) then
                print *, "  Next boundary in direction ", icase, " at distance dc=", dc, " rdir=", rdir
            end if

            isign = rdir_sign(icase)

            rdloc(0:2) = dc(icase) * rdir(0:2) / rdir(icase)
            rdist = sqrt(sum(rdloc**2))

            if (ptau <= ksca_grid * rdist .and. ksca_grid > 0.0_dp) then
                ! scattering event inside current cell
                rdist = ptau / ksca_grid
                rloc = rloc + rdir * rdist

                weight_absorbed = weight * (1.0_dp - exp(-kabs_grid * rdist))
                if (source <= 1) then
                    outradconv(rind(0), rind(1), rind(2), 0) = outradconv(rind(0), rind(1), rind(2), 0) + weight_absorbed
                    outradconv(rind(0), rind(1), rind(2), 1) = outradconv(rind(0), rind(1), rind(2), 1) + weight_absorbed**2
                end if
                
                new_weight = weight - weight_absorbed

                if (source == 2) then
                    phradval = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
                    notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
                    if (debug == 1) then
                        print *, "  phradval", rindsrc(0), rindsrc(1), rindsrc(2), phradval
                    end if
                    outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 0) = outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 0) + phradval*notindsrc
                    outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 1) = outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 1) + (phradval*notindsrc)**2
                else if (source == 4) then
                    if (swlw == 0) then
                        phradval = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
                        outradimg(rindsrc(0), rindsrc(1), 0) = outradimg(rindsrc(0), rindsrc(1), 0) + phradval
                        outradimg(rindsrc(0), rindsrc(1), 1) = outradimg(rindsrc(0), rindsrc(1), 1) + phradval**2
                    else if (swlw == 1) then
                        phrad1val = 1.0_dp
                        r1loc(0:2) = rloc(0:2)
                        r1dir(0:2) = -dirsol(0:2)
                        r1dir_sign(0:2) = max(0, min(1, ceiling(r1dir(0:2))))
                        r1ind(0:2) = rind(0:2)
                        it1 = 0
                        do while (it1 < it1max)
                            if (debug > 0) then
                                if (debug == 1) then
                                    print *, "Local estimation step ", it1, " R1loc: ", r1loc, " R1ind: ", r1ind, " R1dir: ", r1dir
                                else
                                    write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'b', it1, r1loc(0), r1loc(1), r1loc(2), r1ind(0), r1ind(1), r1ind(2), &
                                        r1dir(0), r1dir(1), r1dir(2), phrad1val, 0.0_dp
                                end if
                            end if
                            kext_grid = kext(r1ind(0), r1ind(1), r1ind(2))

                            xbnd(0) = xarr(r1ind(0))
                            xbnd(1) = xarr(r1ind(0)+1)
                            ybnd(0) = yarr(r1ind(1))
                            ybnd(1) = yarr(r1ind(1)+1)
                            zbnd(0) = zarr(r1ind(2))
                            zbnd(1) = zarr(r1ind(2)+1)


                            dc(0) = xbnd(r1dir_sign(0)) - r1loc(0)
                            dc(1) = ybnd(r1dir_sign(1)) - r1loc(1)
                            dc(2) = zbnd(r1dir_sign(2)) - r1loc(2)

                            icase = 2
                            icase = merge(1, icase, abs(dc(1) * r1dir(icase)) < abs(dc(icase) * r1dir(1)))
                            icase = merge(0, icase, abs(dc(0) * r1dir(icase)) < abs(dc(icase) * r1dir(0)))

                            isign = r1dir_sign(icase)

                            rdloc(0:2) = dc(icase) * r1dir(0:2) / r1dir(icase)
                            rdist = sqrt(sum(rdloc**2))

                            ! move to next boundary
                            r1loc(0:2) = r1loc(0:2) + r1dir(0:2) * rdist

                            phrad1val = phrad1val * exp(-kext_grid * rdist)

                            ! update cell index
                            ! r1ind(icase) = r1ind(icase) + isign

                            ! step cell index in icase direction
                            r1ind2(0:2) = r1ind(0:2)
                            r1ind2(icase) = r1ind(icase) + 2 * isign - 1
                            ! wrap x,y
                            r1ind(0) = mod(r1ind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * r1ind(0)
                            r1ind(1) = mod(r1ind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * r1ind(1)
                            r1ind(2) = r1ind2(2)

                            r1loc(icase) = merge( &
                                real(r1ind(icase), dp) * dxs(icase) * real(isign, dp) + real(r1ind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
                                r1loc(icase), &
                                icase <= 1 )
                            it1 = it1 + 1
                            if (r1ind(2) >= nz .or. phrad1val < 1.0e-4_dp) exit
                        end do
                        mu = -sum(dirsol(0:2) * rdir(0:2))
                        g_grid = gparam(rind(0), rind(1), rind(2))
                        phradval = phrad1val * new_weight * (1.0_dp - g_grid**2) &
                                        / (4.0_dp * pi * (1.0_dp + g_grid**2 - 2.0_dp * g_grid * mu)**1.5_dp * abs(dirsol(2)))
                        if (debug == 1) then
                            print *, "  phradval", rindsrc(0), rindsrc(1), rindsrc(2), phradval
                        end if
                        outradimg(rindsrc(0), rindsrc(1), 0) = outradimg(rindsrc(0), rindsrc(1), 0) + phradval
                        outradimg(rindsrc(0), rindsrc(1), 1) = outradimg(rindsrc(0), rindsrc(1), 1) + phradval**2
                    end if

                end if

                if (new_weight < weight_min) then
                    ! Russian Roulette
                    if (debug == 1) then
                        print *, "  Russian Roulette triggered. Weight before RR: ", new_weight
                    end if
                    call random_number(rnd)
                    if (rnd * weight_rr > new_weight) then
                        ! photon terminated
                        if (debug == 1) then
                            print *, "  Photon terminated by RR."
                        end if
                        exit
                    else
                        new_weight = weight_rr
                        if (debug == 1) then
                            print *, "  Photon survived RR. New Weight: ", new_weight
                        end if
                    end if
                end if

                ! Henyey-Greenstein
                g = max(g_grid, 1.0e-8_dp)
                call random_number(rnd); s = 2.0_dp * rnd - 1.0_dp
                mu = (1.0_dp + g**2.0_dp - ((1.0_dp - g**2.0_dp) / (1.0_dp + g * s))**2.0_dp) / (2.0_dp * g)
                call random_number(rnd); phi = 2.0_dp * pi * rnd
                sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
                ux = rdir(0); uy = rdir(1); uz = rdir(2)
                denom = sqrt(ux*ux + uy*uy)
                rdir(0) = merge( &
                    sint * (-ux * uz * cos(phi) + uy * sin(phi)) / denom + ux * mu, &
                    merge(sint * cos(phi), -sint * cos(phi), uz > 0.0_dp), &
                    abs(uz) < 0.9999_dp)

                rdir(1) = merge( &
                    sint * (-uy * uz * cos(phi) - ux * sin(phi)) / denom + uy * mu, &
                    merge(sint * sin(phi), -sint * sin(phi), uz > 0.0_dp), &
                    abs(uz) < 0.9999_dp)

                rdir(2) = merge( &
                    sint * cos(phi) * denom + uz * mu, &
                    merge(mu, -mu, uz > 0.0_dp), &
                    abs(uz) < 0.9999_dp)
                ! normalize
                rdir = rdir / sqrt(sum(rdir**2))
                rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))
                call random_number(rnd); ptau = -log(max(1.0e-300_dp, rnd))
                weight = new_weight
                scaord = scaord + 1
                if (debug == 1) then
                    print *, "  Scattering event. New Rdir: ", rdir, " New Weight: ", weight, 'mu=', mu, 'phi=', phi, 'Ptau=', ptau
                end if
            else
                ! move to next boundary intersection (or remaining segment after scattering)
                ! rloc = rloc + rdir * rdist
                rloc(0:2) = rloc(0:2) + rdir(0:2) * rdist
                ! rloc(0) = mod(rloc(0) + maxx, maxx) * real(isign) + (maxx - mod(maxx - rloc(0), maxx)) * real(1 - isign)
                ! rloc(1) = mod(rloc(1) + maxy, maxy) * real(isign) + (maxy - mod(maxy - rloc(1), maxy)) * real(1 - isign)

                ptau = ptau - ksca_grid * rdist

                weight_absorbed = weight * (1.0_dp - exp(-kabs_grid * rdist))
                if (source <= 1) then
                    outradconv(rind(0), rind(1), rind(2), 0) = outradconv(rind(0), rind(1), rind(2), 0) + weight_absorbed
                    outradconv(rind(0), rind(1), rind(2), 1) = outradconv(rind(0), rind(1), rind(2), 1) + weight_absorbed**2
                else if (source == 2) then
                    phradval = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
                    notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
                    if (debug == 1) then
                        print *, "  phradval", rindsrc(0), rindsrc(1), rindsrc(2), phradval
                    end if
                    outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 0) = outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 0) + phradval*notindsrc
                    outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 1) = outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 1) + (phradval*notindsrc)**2
                end if

                new_weight = weight - weight_absorbed

                ! irradiance sampled along path: new_weight * abs(rdir[icase])
                outradirr(rind(0), rind(1), rind(2), min(min(scaord, 1), swlw), 2 * icase + isign, 0) = &
                    outradirr(rind(0), rind(1), rind(2), min(min(scaord, 1), swlw), 2 * icase + isign, 0) + new_weight * abs(rdir(icase))
                outradirr(rind(0), rind(1), rind(2), min(min(scaord, 1), swlw), 2 * icase + isign, 1) = &
                    outradirr(rind(0), rind(1), rind(2), min(min(scaord, 1), swlw), 2 * icase + isign, 1) + (new_weight * abs(rdir(icase))) ** 2.0_dp

                weight = new_weight

                ! step cell index in icase direction
                rind2(0:2) = rind(0:2)
                rind2(icase) = rind(icase) + 2 * isign - 1
                ! wrap x,y
                rind(0) = mod(rind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * rind(0)
                rind(1) = mod(rind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * rind(1)
                rind(2) = rind2(2)

                if (debug == 1) then
                    print *, "  test", rloc(:)
                end if
                ! rloc(0) = mod(rloc(0) + maxx, maxx) * real(isign) + (maxx - mod(maxx - rloc(0), maxx)) * real(1 - isign)
                ! rloc(1) = mod(rloc(1) + maxy, maxy) * real(isign) + (maxy - mod(maxy - rloc(1), maxy)) * real(1 - isign)
                rloc(icase) = merge( &
                    real(rind(icase), dp) * dxs(icase) * real(isign, dp) + real(rind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
                    rloc(icase), &
                    icase <= 1 )

                if (debug == 1) then
                    print *, "  Moved to next cell. Rloc: ", rloc, " Rind: ", rind, "icase", icase, " Weight: ", weight, 'abs(rdir(icase))=', abs(rdir(icase))
                end if

                it = it + 1

                if (rind(2) >= nz) then
                    if (debug > 0) then
                        if (debug == 1) then
                            print *, "Photon exited TOM."
                        else 
                            write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'a', 0, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), &
                                rdir(0), rdir(1), rdir(2), weight, ptau
                        end if
                    end if
                    exit
                else if (rind(2) < 0) then
                    if (debug > 0) then
                        if (debug == 1) then
                            print *, "Photon hit the ground."
                        else 
                            write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'a', 0, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), &
                                rdir(0), rdir(1), rdir(2), weight, ptau
                        end if
                    end if
                    ! Reflection

                    rind(2) = 0
                    rloc(2) = 0.0_dp + 1.0e-8_dp

                    ! phradval = phradval + weight * (1.0_dp - galb(rind(0), rind(1))) * bgrnd(rind(0), rind(1))
                    new_weight = weight * galb(rind(0), rind(1))

                    if (source == 2) then
                        phradval = weight * (1.0_dp - galb(rind(0), rind(1))) * bgrnd(rind(0), rind(1))
                        if (debug == 1) then
                            print *, "  phradval", rindsrc(0), rindsrc(1), rindsrc(2), phradval
                        end if
                        outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 0) = outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 0) + phradval
                        outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 1) = outradconv(rindsrc(0), rindsrc(1), rindsrc(2), 1) + phradval**2
                    else if (source == 4) then
                        if (swlw == 0) then
                            phradval = weight * (1.0_dp - galb(rind(0), rind(1))) * bgrnd(rind(0), rind(1))
                            outradimg(rindsrc(0), rindsrc(1), 0) = outradimg(rindsrc(0), rindsrc(1), 0) + phradval
                            outradimg(rindsrc(0), rindsrc(1), 1) = outradimg(rindsrc(0), rindsrc(1), 1) + phradval**2
                        else if (swlw == 1) then
                            phrad1val = 1.0_dp
                            r1loc(0:2) = rloc(0:2)
                            r1dir(0:2) = -dirsol(0:2)
                            r1dir_sign(0:2) = max(0, min(1, ceiling(r1dir(0:2))))
                            r1ind(0:2) = rind(0:2)
                            it1 = 0
                            do while (it1 < it1max)
                                if (debug > 0) then
                                    if (debug == 1) then
                                        print *, "Local estimation step ", it1, " R1loc: ", r1loc, " R1ind: ", r1ind, " R1dir: ", r1dir
                                    else
                                        write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'b', it1, r1loc(0), r1loc(1), r1loc(2), r1ind(0), r1ind(1), r1ind(2), &
                                            r1dir(0), r1dir(1), r1dir(2), phrad1val, 0.0_dp
                                    end if
                                end if
                                kext_grid = kext(r1ind(0), r1ind(1), r1ind(2))

                                xbnd(0) = xarr(r1ind(0))
                                xbnd(1) = xarr(r1ind(0)+1)
                                ybnd(0) = yarr(r1ind(1))
                                ybnd(1) = yarr(r1ind(1)+1)
                                zbnd(0) = zarr(r1ind(2))
                                zbnd(1) = zarr(r1ind(2)+1)

                                dc(0) = xbnd(r1dir_sign(0)) - r1loc(0)
                                dc(1) = ybnd(r1dir_sign(1)) - r1loc(1)
                                dc(2) = zbnd(r1dir_sign(2)) - r1loc(2)

                                icase = 2
                                icase = merge(1, icase, abs(dc(1) * r1dir(icase)) < abs(dc(icase) * r1dir(1)))
                                icase = merge(0, icase, abs(dc(0) * r1dir(icase)) < abs(dc(icase) * r1dir(0)))

                                isign = r1dir_sign(icase)

                                rdloc(0:2) = dc(icase) * r1dir(0:2) / r1dir(icase)
                                rdist = sqrt(sum(rdloc**2))

                                ! move to next boundary
                                r1loc(0:2) = r1loc(0:2) + r1dir(0:2) * rdist

                                phrad1val = phrad1val * exp(-kext_grid * rdist)

                                ! update cell index
                                ! r1ind(icase) = r1ind(icase) + isign

                                ! step cell index in icase direction
                                r1ind2(0:2) = r1ind(0:2)
                                r1ind2(icase) = r1ind(icase) + 2 * isign - 1
                                ! wrap x,y
                                r1ind(0) = mod(r1ind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * r1ind(0)
                                r1ind(1) = mod(r1ind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * r1ind(1)
                                r1ind(2) = r1ind2(2)

                                r1loc(icase) = merge( &
                                    real(r1ind(icase), dp) * dxs(icase) * real(isign, dp) + real(r1ind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
                                    r1loc(icase), &
                                    icase <= 1 )

                                it1 = it1 + 1
                                if (r1ind(2) >= nz .or. phrad1val < 1.0e-4_dp) exit
                            end do
                            phradval = phrad1val * new_weight / pi
                            if (debug == 1) then
                                print *, "phradval", phradval
                            end if
                            outradimg(rindsrc(0), rindsrc(1), 0) = outradimg(rindsrc(0), rindsrc(1), 0) + phradval
                            outradimg(rindsrc(0), rindsrc(1), 1) = outradimg(rindsrc(0), rindsrc(1), 1) + phradval**2
                        end if

                    end if

                    if (new_weight < weight_min) then
                        ! Russian Roulette
                        if (debug == 1) then
                            print *, "  Russian Roulette triggered. Weight before RR: ", new_weight
                        end if
                        call random_number(rnd)
                        if (rnd * weight_rr > new_weight) then
                            ! photon terminated
                            if (debug == 1) then
                                print *, "  Photon terminated by RR."
                            end if
                            exit
                        else
                            new_weight = weight_rr
                            if (debug == 1) then
                                print *, "  Photon survived RR. New Weight: ", new_weight
                            end if
                        end if
                    end if
                    
                    call random_number(rnd); mu = sqrt(rnd)
                    call random_number(rnd); phi = 2.0_dp * pi * rnd
                    sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
                    rdir(0) = sint * cos(phi)
                    rdir(1) = sint * sin(phi)
                    rdir(2) = mu
                    ! normalize
                    rdir = rdir / sqrt(sum(rdir**2))
                    rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))

                    call random_number(rnd); ptau = -log(max(1.0e-300_dp, rnd))
                    scaord = scaord + 1
                    weight = new_weight
                end if
            end if
        end do
    end do
    end do
    end do
    end do
    end do

    if (debug >= 2) then
        close(iutraj)
    end if

    write(*,*) "All photons processed. Writing..."

    do ix = 0, nx-1
        do iy = 0, ny-1
            do iz = 0, nz-1
                do idi = 0, swlw
                    do icase = 0, 5
                        outradirr(ix,iy,iz,idi,icase,0) = outradirr(ix,iy,iz,idi,icase,0) / real(nphoton, dp)
                        outradirr(ix,iy,iz,idi,icase,1) = outradirr(ix,iy,iz,idi,icase,1) / real(nphoton, dp)
                    end do
                end do
                outradconv(ix,iy,iz,0) = (outradconv(ix,iy,iz,0) / real(nphoton, dp) &
                    - real(1 - swlw) * 1.0_dp * bplnk(ix,iy,iz)) * kabs(ix,iy,iz)
                    ! ) * kabs(ix,iy,iz)
                outradconv(ix,iy,iz,1) = (outradconv(ix,iy,iz,1) / real(nphoton, dp)) &
                    * ((real(1 - swlw) * 1.0_dp * bplnk(ix,iy,iz)) * kabs(ix,iy,iz))**2.0_dp
                    ! * kabs(ix,iy,iz)
            end do
        outradimg(ix,iy,0) = outradimg(ix,iy,0) / real(nphoton, dp)
        outradimg(ix,iy,1) = outradimg(ix,iy,1) / real(nphoton, dp)
        end do
    end do

    call write_output_d6(outradirr(0:nx-1,0:ny-1,0:nz-1,0:swlw,4:5,0:1), nx, ny, nz, swlw+1, 2, nphoton, trim(line)//"/outradirr.txt")
    call write_output_d4(outradconv(0:nx-1,0:ny-1,0:nz-1,0:1), nx, ny, nz, nphoton, trim(line)//"/outradconv.txt")
    call write_output_d3(outradimg(0:nx-1,0:ny-1,0:1), nx, ny, nphoton, trim(line)//"/outradimg.txt")

    write(*,*) "Done."
    
contains

    subroutine write_output_d6(outrad, nx, ny, nz, nd, na, nphoton, filename)
        implicit none
        real(dp), intent(in) :: outrad(0:nx-1,0:ny-1,0:nz-1,0:nd-1,0:na-1,0:1)
        integer, intent(in) :: nx, ny, nz, nd, na, nphoton
        character(len=*), intent(in) :: filename
        real(dp) :: val, val2
        integer :: unit, ix, iy, iz, idi, ia

        unit = 10
        open(unit, file = filename, status = "replace", action = "write", form = "formatted")

        write(unit, '(A)') "! nx ny nz nd na comp nphoton"
        write(unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I10)') nx, ny, nz, nd, na, 2, nphoton
        do ix = 0, nx-1
            do iy = 0, ny-1
                do iz = 0, nz-1
                    do idi = 0, nd-1
                        do ia = 0, na-1
                            val = outrad(ix,iy,iz,idi,ia,0)
                            val2 = (outrad(ix,iy,iz,idi,ia,1) - val**2.0_dp) / real(nphoton, dp)
                            if (abs(val) < 1.0e-20_dp) val = 0.0_dp
                            if (abs(val2) < 1.0e-20_dp) val2 = 0.0_dp
                            write(unit, '(E15.6, E15.6)') val, val2
                        end do
                    end do
                end do
            end do
        end do

        close(unit)
    end subroutine write_output_d6

    subroutine write_output_d4(outrad, nx, ny, nz, nphoton, filename)
        implicit none
        real(dp), intent(in) :: outrad(0:nx-1,0:ny-1,0:nz-1,0:1)
        integer, intent(in) :: nx, ny, nz, nphoton
        character(len=*), intent(in) :: filename
        real(dp) :: val, val2
        integer :: unit, ix, iy, iz

        unit = 10
        open(unit, file = filename, status = "replace", action = "write", form = "formatted")

        ! write(unit, '(A)') "! nx ny nz ncase comp nphoton"
        ! write(unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I10)') nx, ny, nz, 3, 2, nphoton
        write(unit, '(A)') "! nx ny nz comp nphoton"
        write(unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I10)') nx, ny, nz, 2, nphoton
        do ix = 0, nx-1
            do iy = 0, ny-1
                do iz = 0, nz-1
                    ! do icase = 0, 2
                        ! val = outrad(ix,iy,iz,icase,0)
                        ! val2 = (outrad(ix,iy,iz,icase,1) - val**2.0_dp) / real(nphoton, dp)
                        val = outrad(ix,iy,iz,0)
                        val2 = (outrad(ix,iy,iz,1) - val**2.0_dp) / real(nphoton, dp)
                        if (abs(val) < 1.0e-20_dp) val = 0.0_dp
                        if (abs(val2) < 1.0e-20_dp) val2 = 0.0_dp
                        write(unit, '(E15.6, E15.6)') val, val2
                    ! end do
                end do
            end do
        end do

        close(unit)
    end subroutine write_output_d4

    subroutine write_output_d3(outrad, nx, ny, nphoton, filename)
        implicit none
        real(dp), intent(in) :: outrad(0:nx-1,0:ny-1,0:1)
        integer, intent(in) :: nx, ny, nphoton
        character(len=*), intent(in) :: filename
        real(dp) :: val, val2
        integer :: unit, ix, iy

        unit = 10
        open(unit, file = filename, status = "replace", action = "write", form = "formatted")

        write(unit, '(A)') "! nx ny comp nphoton"
        write(unit, '(I5,1X,I5,1X,I5,1X,I10)') nx, ny, 2, nphoton
        do ix = 0, nx-1
            do iy = 0, ny-1
                val = outrad(ix,iy,0)
                val2 = outrad(ix,iy,1)
                if (abs(val) < 1.0e-20_dp) val = 0.0_dp
                if (abs(val2) < 1.0e-20_dp) val2 = 0.0_dp
                write(unit, '(E15.6, E15.6)') val, val2
            end do
        end do

        close(unit)
    end subroutine write_output_d3

end program test_fullmc




