program test_fullmc
    use fullmc_funcs
    implicit none
    ! integer, parameter :: dp = selected_real_kind(15, 307)
    character(len=256) :: line
    integer  :: nphoton
    integer  :: nx, ny, nz
    real(dp) :: dx, dy, dz
    real(dp) :: dxs(0:2)
    real(dp), allocatable :: xarr(:), yarr(:), zarr(:)
    real(dp), allocatable :: kext(:,:,:), ksca(:,:,:), kabs(:,:,:), gparam(:,:,:)
    real(dp), allocatable :: galb(:,:), bplnk(:,:,:), bgrnd(:,:)
    integer  :: source, swlw, transfer_mode
    real(dp) :: solmu, solphi
    real(dp) :: viewmu, viewphi
    integer  :: seedval
    integer  :: scaord
    real(dp) :: dirsol(0:2), dirview(0:2)
    integer  :: iphoton, it, itmax, itlpbmax
    ! integer  :: itlpb
    real(dp) :: maxx, maxy, maxz
    real(dp) :: rloc(0:2), rdir(0:2)
    integer  :: rind(0:2), rindsrc(0:2)
    ! real(dp) :: rlpbloc(0:2), rlpbdir(0:2)
    ! integer  :: rlpbind(0:2)
    ! integer  :: rind2(0:2), rlpbind2(0:2)
    integer  :: rdir_sign(0:2)
    ! integer  :: rlpbdir_sign(0:2)
    ! real(dp) :: xbnd(0:1), ybnd(0:1), zbnd(0:1), dc(0:2)
    real(dp) :: rdloc(0:2)
    real(dp) :: init_weight
    real(dp) :: weight, ptau, rnd, vqla, vqllpb
    real(dp) :: ksca_grid, kabs_grid, g_grid
    ! real(dp) :: kext_grid
    integer  :: ix, iy, iz, ia
    integer  :: nzmin, nzmax, namax
    integer  :: icase, isign, idi
    real(dp) :: notindsrc
    real(dp) :: rdist
    real(dp) :: weight_absorbed, new_weight
    real(dp) :: weight_min, weight_rr
    logical  :: survived
    real(dp) :: g, mu, phi, sint
    ! real(dp) :: s, ux, uy, uz, denom
    integer  :: recind(0:2), irx, iry, irz, irdi, ira
    real(dp) :: nphotot
    real(dp) :: escale
    procedure(rec_scat_iface), pointer :: recorder_scattering => null()
    procedure(rec_bound_iface), pointer :: recorder_boundary => null()

    real(dp), allocatable :: recflx(:,:,:,:,:,:), recconv(:,:,:,:), recimg(:,:,:)
    real(dp), allocatable :: phflx(:,:,:,:,:), phconv(:,:,:), phimg(:,:)
    real(dp), allocatable :: outflx(:,:,:,:,:,:), outconv(:,:,:,:), outimg(:,:,:)
    real(dp), allocatable :: outflx_tmp(:,:,:,:,:,:)
    real(dp) :: flx0, flx1, conv0, conv1, img0, img1
    ! real(dp) :: pi = acos(-1.0_dp)
    integer  :: i
    integer  :: iuconf = 41
    integer  :: iutraj = 42
    integer  :: seed_size
    integer, allocatable :: seed(:)
    integer  :: debug
    character(len=256) :: dbgmsg

    !-- Read configuration file name
    if (command_argument_count() >= 1) then
        call get_command_argument(1, line)
    else
        read(*,'(A)') line
    end if

    !-- Read configuration file contents
    ! open(unit=iuconf, file='out/test_fullmc/config.txt', status='old', action='read')
    open(unit=iuconf, file=trim(line)//'/config.txt', status='old', action='read')
    read(iuconf,*) nx, ny, nz
    read(iuconf,*) dx, dy, dz
    read(iuconf,*) source, swlw
    read(iuconf,*) transfer_mode
    read(iuconf,*) solmu, solphi
    read(iuconf,*) viewmu, viewphi
    read(iuconf,*) nphoton
    read(iuconf,*) seedval
    read(iuconf,*) debug

    write(*,*) 'Configuration: nx=', nx, ' ny=', ny, ' nz=', nz
    write(*,*) '               dx=', dx, ' dy=', dy, ' dz=', dz
    write(*,*) '               nphoton=', nphoton
    write(*,*) '               source=', source
    write(*,*) '               swlw=', swlw
    write(*,*) '               transfer_mode=', transfer_mode
    write(*,*) '               seedval=', seedval
    write(*,*) '               debug=', debug

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
            end do
            read(iuconf,*) galb(ix,iy), bgrnd(ix,iy)
        end do
    end do
    !-- End of configuration file reading

    allocate(xarr(0:nx), yarr(0:ny), zarr(0:nz))
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

    init_weight = 1.0_dp
    weight_min = 0.5_dp
    weight_rr = 0.5_dp

    dirsol(0) = sqrt(1.0_dp - solmu**2) * sin(solphi)
    dirsol(1) = sqrt(1.0_dp - solmu**2) * cos(solphi)
    dirsol(2) = -solmu
    dirsol = dirsol / sqrt(sum(dirsol**2))

    dirview(0) = sqrt(1.0_dp - viewmu**2) * sin(viewphi)
    dirview(1) = sqrt(1.0_dp - viewmu**2) * cos(viewphi)
    dirview(2) = -viewmu
    dirview = dirview / sqrt(sum(dirview**2))

    ! RNG seeding for reproducibility
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = seedval + (/ (i-1, i=1,seed_size) /)
    call random_seed(put=seed)

    maxx = real(nx,dp) * dx
    maxy = real(ny,dp) * dy
    maxz = real(nz,dp) * dz
    dxs(0) = dx
    dxs(1) = dy
    dxs(2) = dz

    if (source <= 1) then
        recorder_scattering => recorder_scat1
        recorder_boundary => recorder_bound1
    else if (source == 2) then
        recorder_scattering => recorder_scat2
        recorder_boundary => recorder_bound2
    else if (source == 3) then
        recorder_scattering => recorder_scat3
        recorder_boundary => recorder_bound3
    else if (source == 4) then
        recorder_scattering => recorder_scat4
        recorder_boundary => recorder_bound4
    end if

    allocate(recflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5,0:1))
    allocate(recconv(0:nx-1,0:ny-1,0:nz-1,0:1))
    allocate(recimg(0:nx-1,0:ny-1,0:1))
    recflx(:, :, :, :, :, :) = 0.0_dp
    recconv(:, :, :, :) = 0.0_dp
    recimg(:, :, :) = 0.0_dp

    allocate(phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5))
    allocate(phconv(0:nx-1,0:ny-1,0:nz-1))
    allocate(phimg(0:nx-1,0:ny-1))

    itmax = 100000
    itlpbmax = 100000

    if (debug >= 2) then
        open (unit=iutraj, file=trim(line)//'/photon_trajectory.txt', status='replace', action='write')
        write(iutraj,*) 'ix, iy, iz, ia, iphoton, it, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), '// &
                    'rdir(0), rdir(1), rdir(2), weight, ptau'
    end if

    write(*,*) "Starting photon transport..."

    ! nzmax = merge(merge(0, nz - 1, source <= 1), 0, source <= 3)
    ! namax = merge(0, 5, source /= 3)
    if (source <= 1) then
        ! Insert from TOA
        nzmin = 0
        nzmax = 0
        namax = 0
    else if (source == 2) then
        ! Insert from within grid volume
        nzmin = 0
        nzmax = nz - 1
        namax = 0
    else if (source == 3) then
        ! Insert from grid sides
        nzmin = -1
        nzmax = nz - 1
        namax = 5
    else if (source == 4) then
        ! Insert from TOA for satellite imaging
        nzmin = 0
        nzmax = 0
        namax = 0
    end if

    do ix = 0, nx - 1 ! X
    do iy = 0, ny - 1 ! Y
    do iz = nzmin, nzmax  ! Z
    do ia = 0, namax  ! Sides

    if ((source == 3) .and. (iz < 0) .and. (ia /= 5)) cycle ! top sources only at ground level

    do iphoton = 0, nphoton - 1
    
    !-- Photon initialization
    !-- location
    call photon_initloc(source, ix, iy, iz, ia, xarr, yarr, zarr, dx, dy, dz, nx, ny, nz, maxx, maxy, maxz, &
        rloc, rindsrc, rind)
    !-- direction
    call photon_initdir(source, dirsol, dirview, rdir, rdir_sign, ia)
    !-- weight
    weight = init_weight
    !-- optical depth to next event
    call random_number(rnd); ptau = -log(max(1.0e-300_dp, rnd))
    !-- scattering order
    scaord = 0

    phflx(:, :, :, :, :) = 0.0_dp
    phconv(:, :, :) = 0.0_dp
    phimg(:, :) = 0.0_dp

    !<< Sampling >>
    if (source <= 1) then
        vqla = weight !* abs(rdir(2))
        call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), 4, nx, ny, nz, swlw)
    end if

    !-- Main photon transport iteration (it: number of grid crossing)
    it = 0
    do while (it < itmax)

        if (debug == 1) write(*,*) "Step ", it, " Rloc: ", rloc, " Rind: ", rind, " Rdir: ", rdir, " Weight: ", weight, " Ptau: ", ptau
        if (debug == 2) call record_trajectory(iutraj, ix, iy, iz, ia, iphoton, it, 'a', 0, rloc, rind, rdir, weight, ptau)

        ksca_grid = ksca(rind(0), rind(1), rind(2))
        kabs_grid = kabs(rind(0), rind(1), rind(2))
        g_grid = gparam(rind(0), rind(1), rind(2))

        call photon_intersect(nx, ny, nz, rind, rloc, rdir, rdir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)

        if (ptau <= ksca_grid * rdist .and. ksca_grid > 0.0_dp) then ! scattering event inside current cell
            
            rdist = ptau / ksca_grid
            rloc = rloc + rdir * rdist
            weight_absorbed = weight * (1.0_dp - exp(-kabs_grid * rdist))
            new_weight = weight - weight_absorbed

            !<< Sampling >>
            ! call sample_scattering(weight_absorbed, new_weight, source, ia, scaord, swlw, &
            !     nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
            !     rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
            !     ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
            call recorder_scattering(weight_absorbed, new_weight, ia, scaord, swlw, &
                nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
                rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
                ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
            ! if (source <= 1) then
            !     vqla = weight_absorbed
            !     call sample_d4(vqla, phconv, rind, nx, ny, nz)
            ! else if (source == 2) then
            !     notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
            !     vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * notindsrc
            !     call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)
            ! else if (source == 3) then
            !     vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * 0.5_dp
            !     call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)
            ! else if (source == 4) then
            !     if (swlw == 0) then
            !         vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
            !         ! call sample_d3(vqla, phimg, rindsrc, nx, ny)
            !     else if (swlw == 1) then
            !         call photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
            !             ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
            !         mu = -sum(dirsol(0:2) * rdir(0:2))
            !         g_grid = gparam(rind(0), rind(1), rind(2))
            !         vqla = vqllpb * new_weight * (1.0_dp - g_grid**2) &
            !             / (4.0_dp * pi * (1.0_dp + g_grid**2 - 2.0_dp * g_grid * mu)**1.5_dp * abs(dirsol(2)))
            !     end if
            !     call sample_d3(vqla, phimg, rindsrc, nx, ny)
            ! end if

            call photon_rroulette(new_weight, weight_min, weight_rr, survived, debug)
            if (.not. survived) exit
            weight = new_weight
            call photon_scattering(g_grid, rdir, rdir_sign, ptau, scaord)

        else ! move to next boundary intersection

            ptau = ptau - ksca_grid * rdist
            weight_absorbed = weight * (1.0_dp - exp(-kabs_grid * rdist))
            new_weight = weight - weight_absorbed

            !<< Sampling >>
            call recorder_boundary(weight_absorbed, new_weight, ia, scaord, swlw, icase, isign, rdir, &
            rind, rindsrc, nx, ny, nz, phflx, phconv, bplnk)
            ! if (source <= 1) then
            !     vqla = new_weight !* abs(rdir(icase))
            !     call sample_d6(vqla, phflx, rind, min(min(scaord, 1), swlw), 2 * icase + isign, nx, ny, nz, swlw)
            !     vqla = weight_absorbed
            !     call sample_d4(vqla, phconv, rind, nx, ny, nz)
            ! else if (source == 2) then
            !     notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
            !     vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * notindsrc
            !     call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)
            ! else if (source == 3) then
            !     vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * 0.5_dp
            !     call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)
            ! end if

            call photon_rroulette(new_weight, weight_min, weight_rr, survived, debug)
            if (.not. survived) exit

            weight = new_weight
            call photon_movegrid(rind, rloc, rdir, rdist, icase, isign, nx, ny, dxs, transfer_mode)
            it = it + 1

            if (rind(2) >= nz) then ! TOA
                if (debug == 1) write(*,*) "Photon exited TOM."
                if (debug == 2) call record_trajectory(iutraj, ix, iy, iz, ia, iphoton, it, 'a', 0, rloc, rind, rdir, weight, ptau)
                exit ! TOA Exit

            else if (rind(2) < 0) then ! Ground
                if (debug == 1) write(*,*) "Photon hit the ground."
                if (debug == 2) call record_trajectory(iutraj, ix, iy, iz, ia, iphoton, it, 'a', 0, rloc, rind, rdir, weight, ptau)

                new_weight = weight * galb(rind(0), rind(1))

                !<< Sampling >>
                if (source == 2) then
                    vqla = weight * (1.0_dp - galb(rind(0), rind(1))) * bgrnd(rind(0), rind(1))
                    call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)
                else if (source == 3) then
                    vqla = weight * (1.0_dp - galb(rind(0), rind(1))) * bgrnd(rind(0), rind(1)) * 0.5_dp
                    call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)
                else if (source == 4) then
                    if (swlw == 0) then
                        vqla = weight * (1.0_dp - galb(rind(0), rind(1))) * bgrnd(rind(0), rind(1))
                        ! call sample_d3(vqla, phimg, rindsrc, nx, ny)
                    else if (swlw == 1) then
                        rind(2) = 0
                        rloc(2) = 0.0_dp
                        ! rloc(2) = 0.0_dp + 1.0e-8_dp
                        call photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
                            ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
                        vqla = vqllpb * new_weight / pi
                    end if
                    call sample_d3(vqla, phimg, rindsrc, nx, ny)
                end if

                call photon_rroulette(new_weight, weight_min, weight_rr, survived, debug)
                if (.not. survived) exit
                call photon_reflection(rdir, rdir_sign, ptau, scaord)
                weight = new_weight

                !<< Sampling >>
                if (source <= 1) then
                    ! rind(2) = -1
                    vqla = weight !* abs(rdir(2))
                    call sample_d6(vqla, phflx, rind, min(min(scaord, 1), swlw), 5, nx, ny, nz, swlw)
                end if

                rind(2) = 0
                rloc(2) = 0.0_dp
                ! rloc(2) = 0.0_dp + 1.0e-8_dp

            end if
        end if
    end do ! while (it < itmax)
    if (it >= itmax) then
        write(dbgmsg, '(A,I10)') "Warning: Photon transport reached maximum iteration at photon ", iphoton
        write(*,*) trim(dbgmsg)
    end if
    !-- Photon sample recording
    if (source <= 1) then
        do irx = 0, nx - 1
        do iry = 0, ny - 1
        recind(0) = irx
        recind(1) = iry
        do irz = -1, nz
            recind(2) = irz
            do irdi = 0, swlw
                do ira = 0, 5
                    call store_d6(phflx, recflx, recind, irdi, ira, nx, ny, nz, swlw)
                end do
            end do
        end do ! iz
        do irz = 0, nz - 1
            recind(2) = irz
            call store_d4(phconv, recconv, recind, nx, ny, nz)
        end do ! iz
        end do ! iy
        end do ! ix
    else if (source == 2) then
        recind(0:2) = (/ix, iy, iz/)
        call store_d4(phconv, recconv, recind, nx, ny, nz)
    else if (source == 3) then
        recind(0:2) = (/ix, iy, iz/)
        do irdi = 0, swlw
            call store_d6(phflx, recflx, recind, irdi, ia, nx, ny, nz, swlw)
        end do
    else if (source == 4) then
        recind(0:2) = (/ix, iy, 0/)
        call store_d3(phimg, recimg, recind, nx, ny)
    end if
    end do ! iphoton
    end do ! ia
    end do ! iz
    end do ! iy
    end do ! ix

    write(*,*) "Photon transport completed."

    if (debug >= 2) close(iutraj)

    write(*,*) "All photons processed. Writing..."

    allocate(outflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5,0:1))
    allocate(outconv(0:nx-1,0:ny-1,0:nz-1,0:1))
    allocate(outimg(0:nx-1,0:ny-1,0:1))

    write(*,*) "Normalizing results..."

    if (source <= 1) then
        escale = real(nx * ny, dp)
    else
        escale = 1.0_dp
    end if
    nphotot = real(nphoton, dp) * escale

    do ix = 0, nx-1
    do iy = 0, ny-1
    do iz = -1, nz
        do idi = 0, swlw
            do icase = 0, 5
                ! recflx(ix,iy,iz,idi,icase,0) = recflx(ix,iy,iz,idi,icase,0) / real(nphoton, dp)
                ! recflx(ix,iy,iz,idi,icase,1) = recflx(ix,iy,iz,idi,icase,1) / real(nphoton, dp)
                flx0 = recflx(ix,iy,iz,idi,icase,0) * escale / nphotot
                flx1 = (recflx(ix,iy,iz,idi,icase,1) * escale**2 / nphotot - flx0**2) / (nphotot - 1.0_dp)
                ! flx1 = recflx(ix,iy,iz,idi,icase,1) / real(nphoton, dp)
                outflx(ix,iy,iz,idi,icase,0) = flx0
                outflx(ix,iy,iz,idi,icase,1) = flx1
            end do
        end do
    end do ! iz
    end do ! iy
    end do ! ix

    write(*,*) "Calculating convergence..."
    
    if (source == 3) then
        do ix = 0, nx-1
        do iy = 0, ny-1
        do iz = 0, nz-1
            ! recconv(ix,iy,iz,0:1) = 0.0_dp
            conv0 = 0.0_dp
            conv1 = 0.0_dp
            do idi = 0, swlw
                ! recconv(ix,iy,iz,0) = &
                ! - recflx(ix,iy,iz,idi,0,0) * dy * dz &
                ! - recflx(ix,iy,iz,idi,1,0) * dy * dz &
                ! - recflx(ix,iy,iz,idi,2,0) * dx * dz &
                ! - recflx(ix,iy,iz,idi,3,0) * dx * dz &
                ! - recflx(ix,iy,iz,idi,4,0) * dx * dy &
                ! - recflx(ix,iy,iz,idi,5,0) * dx * dy &
                ! + recflx(mod(ix - 1 + nx, nx)*transfer_mode+ix*(1 - transfer_mode),iy,iz,idi,1,0) * dy * dz &
                ! + recflx(mod(ix + 1 + nx, nx)*transfer_mode+ix*(1 - transfer_mode),iy,iz,idi,0,0) * dy * dz &
                ! + recflx(ix,mod(iy - 1 + ny, ny)*transfer_mode+iy*(1 - transfer_mode),iz,idi,3,0) * dx * dz &
                ! + recflx(ix,mod(iy + 1 + ny, ny)*transfer_mode+iy*(1 - transfer_mode),iz,idi,2,0) * dx * dz &
                ! + recflx(ix,iy,iz - 1,idi,5,0) * dx * dy &
                ! + recflx(ix,iy,iz + 1,idi,4,0) * dx * dy
                conv0 = conv0 + &
                        (- outflx(ix,iy,iz,idi,0,0) * dy * dz &
                         - outflx(ix,iy,iz,idi,1,0) * dy * dz &
                         - outflx(ix,iy,iz,idi,2,0) * dx * dz &
                         - outflx(ix,iy,iz,idi,3,0) * dx * dz &
                         - outflx(ix,iy,iz,idi,4,0) * dx * dy &
                         - outflx(ix,iy,iz,idi,5,0) * dx * dy &
                         + outflx(mod(ix - 1 + nx, nx)*transfer_mode+ix*(1 - transfer_mode),iy,iz,idi,1,0) * dy * dz &
                         + outflx(mod(ix + 1 + nx, nx)*transfer_mode+ix*(1 - transfer_mode),iy,iz,idi,0,0) * dy * dz &
                         + outflx(ix,mod(iy - 1 + ny, ny)*transfer_mode+iy*(1 - transfer_mode),iz,idi,3,0) * dx * dz &
                         + outflx(ix,mod(iy + 1 + ny, ny)*transfer_mode+iy*(1 - transfer_mode),iz,idi,2,0) * dx * dz &
                         + outflx(ix,iy,iz - 1,idi,5,0) * dx * dy &
                         + outflx(ix,iy,iz + 1,idi,4,0) * dx * dy ) &
                         / (2.0_dp * dx * dy * dz)
                conv1 = conv1 + &
                        (  outflx(ix,iy,iz,idi,0,1) * (dy * dz)**2 &
                         + outflx(ix,iy,iz,idi,1,1) * (dy * dz)**2 &
                         + outflx(ix,iy,iz,idi,2,1) * (dx * dz)**2 &
                         + outflx(ix,iy,iz,idi,3,1) * (dx * dz)**2 &
                         + outflx(ix,iy,iz,idi,4,1) * (dx * dy)**2 &
                         + outflx(ix,iy,iz,idi,5,1) * (dx * dy)**2 &
                         + outflx(mod(ix - 1 + nx, nx)*transfer_mode+ix*(1 - transfer_mode),iy,iz,idi,1,1) * (dy * dz)**2 &
                         + outflx(mod(ix + 1 + nx, nx)*transfer_mode+ix*(1 - transfer_mode),iy,iz,idi,0,1) * (dy * dz)**2 &
                         + outflx(ix,mod(iy - 1 + ny, ny)*transfer_mode+iy*(1 - transfer_mode),iz,idi,3,1) * (dx * dz)**2 &
                         + outflx(ix,mod(iy + 1 + ny, ny)*transfer_mode+iy*(1 - transfer_mode),iz,idi,2,1) * (dx * dz)**2 &
                         + outflx(ix,iy,iz - 1,idi,5,1) * (dx * dy)**2 &
                         + outflx(ix,iy,iz + 1,idi,4,1) * (dx * dy)**2 ) &
                         / ( (2.0_dp * dx * dy * dz)**2 )
            end do
            outconv(ix,iy,iz,0) = conv0
            outconv(ix,iy,iz,1) = conv1
        end do ! iz
        end do ! iy
        end do ! ix
    else
        do ix = 0, nx-1
        do iy = 0, ny-1
        do iz = 0, nz-1
            ! recconv(ix,iy,iz,0) = (recconv(ix,iy,iz,0) / real(nphoton, dp) &
            !     - real(1 - swlw) * 1.0_dp * bplnk(ix,iy,iz)) * kabs(ix,iy,iz)
            !     ! ) * kabs(ix,iy,iz)
            ! recconv(ix,iy,iz,1) = (recconv(ix,iy,iz,1) / real(nphoton, dp)) &
            !     * ((real(1 - swlw) * 1.0_dp * bplnk(ix,iy,iz)) * kabs(ix,iy,iz))**2.0_dp
            !     ! * kabs(ix,iy,iz)
            conv0 = recconv(ix,iy,iz,0) * escale / nphotot
            conv1 = (recconv(ix,iy,iz,1) * escale**2 / nphotot - conv0**2) / (nphotot - 1.0_dp)
            outconv(ix,iy,iz,0) = (conv0 - real(1 - swlw) * 1.0_dp * bplnk(ix,iy,iz)) * kabs(ix,iy,iz)
            outconv(ix,iy,iz,1) = conv1 * kabs(ix,iy,iz)**2.0_dp
        end do ! iz
        end do ! iy
        end do ! ix
    end if

    write(*,*) "Calculating images..."

    do ix = 0, nx-1
    do iy = 0, ny-1
    ! recimg(ix,iy,0) = recimg(ix,iy,0) / real(nphoton, dp)
    ! recimg(ix,iy,1) = recimg(ix,iy,1) / real(nphoton, dp)
    img0 = recimg(ix,iy,0) * escale / nphotot
    img1 = (recimg(ix,iy,1) * escale**2 / nphotot - img0**2) / (nphotot - 1.0_dp)
    outimg(ix,iy,0) = img0
    outimg(ix,iy,1) = img1
    end do ! iy
    end do ! ix

    deallocate(kext)
    deallocate(ksca)
    deallocate(kabs)
    deallocate(gparam)
    deallocate(galb)
    deallocate(bplnk)
    deallocate(bgrnd)
    deallocate(xarr)
    deallocate(yarr)
    deallocate(zarr)
    deallocate(seed)

    deallocate(recflx)
    deallocate(recconv)
    deallocate(recimg)

    write(*,*) "Preparing output arrays..."

    allocate(outflx_tmp(0:nx-1,0:ny-1,0:nz-1,0:swlw,0:1,0:1))
    
    do ix = 0, nx-1
    do iy = 0, ny-1
    do iz = 0, nz-1
        do idi = 0, swlw
            outflx_tmp(ix,iy,iz,idi,0,0) = outflx(ix,iy,iz,idi,4,0)
            outflx_tmp(ix,iy,iz,idi,0,1) = outflx(ix,iy,iz,idi,4,1)
            outflx_tmp(ix,iy,iz,idi,1,0) = outflx(ix,iy,iz,idi,5,0)
            outflx_tmp(ix,iy,iz,idi,1,1) = outflx(ix,iy,iz,idi,5,1)
        end do
    end do
    end do
    end do

    deallocate(outflx)

    write(*,*) "Writing output files to ", trim(line)

    ! call write_output_d6(outflx(0:nx-1,0:ny-1,0:nz-1,0:swlw,4:5,0:1), nx, ny, nz, swlw+1, 2, nphoton, trim(line)//"/outflx.txt")
    call write_output_d6(outflx_tmp(0:nx-1,0:ny-1,0:nz-1,0:swlw,0:1,0:1), nx, ny, nz, swlw+1, 2, nphoton, trim(line)//"/outradirr.txt")
    ! call write_output_d6(outflx(0:nx-1,0:ny-1,0:nz-1,0:swlw,0:1,0:1), nx, ny, nz, swlw+1, 2, nphoton, trim(line)//"/outflx.txt")
    call write_output_d4(outconv(0:nx-1,0:ny-1,0:nz-1,0:1), nx, ny, nz, nphoton, trim(line)//"/outradconv.txt")
    call write_output_d3(outimg(0:nx-1,0:ny-1,0:1), nx, ny, nphoton, trim(line)//"/outradimg.txt")

    deallocate(outflx_tmp)
    deallocate(outconv)
    deallocate(outimg)

    write(*,*) "Done."
    
! contains

    ! subroutine photon_initloc(source, ix, iy, iz, ia, xarr, yarr, zarr, dx, dy, dz, nx, ny, nz, maxx, maxy, maxz, rloc, rindsrc, rind)
    !     implicit none
    !     integer, intent(in) :: source, ix, iy, iz, ia, nx, ny, nz
    !     real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
    !     real(dp), intent(in) :: dx, dy, dz
    !     real(dp), intent(in) :: maxx, maxy, maxz
    !     real(dp), intent(out) :: rloc(0:2)
    !     integer, intent(out) :: rindsrc(0:2), rind(0:2)
    !     real(dp) :: rnd
    !     if (source <= 1) then
    !         ! Top of model
    !         call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
    !         call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
    !         rloc(0) = mod(rloc(0) + maxx, maxx)
    !         rloc(1) = mod(rloc(1) + maxy, maxy)
    !         rloc(2) = maxz

    !         rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !         rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !         rindsrc(2) = nz
    !         rind(0) = rindsrc(0)
    !         rind(1) = rindsrc(1)
    !         rind(2) = rindsrc(2) - 1
    !     else if (source == 2) then
    !         ! Volumetric source
    !         call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
    !         call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
    !         call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
    !         rloc(0) = mod(rloc(0) + maxx, maxx)
    !         rloc(1) = mod(rloc(1) + maxy, maxy)
    !         rloc(2) = mod(rloc(2) + maxz, maxz)
            
    !         rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !         rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !         rindsrc(2) = mod(int( rloc(2) / dz ), nz)
    !         rind(0:2) = rindsrc(0:2)
    !     else if (source == 3) then
    !         ! Grid sides
    !         if (ia == 0) then ! x+
    !             rloc(0) = xarr(ix)
    !             call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
    !             call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
    !             rindsrc(0) = mod(int((rloc(0) + 1.0e-8_dp) / dx ), nx)
    !             rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !             rindsrc(2) = mod(int( rloc(2) / dz ), nz)
    !         else if (ia == 1) then ! x-
    !             rloc(0) = xarr(ix + 1)
    !             call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
    !             call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
    !             rindsrc(0) = mod(int((rloc(0) - 1.0e-8_dp) / dx ), nx)
    !             rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !             rindsrc(2) = mod(int( rloc(2) / dz ), nz)
    !         else if (ia == 2) then ! y+
    !             call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
    !             rloc(1) = yarr(iy)
    !             call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
    !             rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !             rindsrc(1) = mod(int((rloc(1) + 1.0e-8_dp) / dy ), ny)
    !             rindsrc(2) = mod(int( rloc(2) / dz ), nz)
    !         else if (ia == 3) then ! y-
    !             call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
    !             rloc(1) = yarr(iy + 1)
    !             call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
    !             rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !             rindsrc(1) = mod(int((rloc(1) - 1.0e-8_dp) / dy ), ny)
    !             rindsrc(2) = mod(int( rloc(2) / dz ), nz)
    !         else if (ia == 4) then ! z+
    !             call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
    !             call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
    !             rloc(2) = zarr(iz)
    !             rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !             rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !             rindsrc(2) = int(floor((rloc(2) + 1.0e-8_dp) / dz ))
    !         else if (ia == 5) then ! z-
    !             call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
    !             call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
    !             rloc(2) = zarr(iz + 1)
    !             rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !             rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !             rindsrc(2) = int(floor((rloc(2) - 1.0e-8_dp) / dz ))
    !         end if
    !         rind(0:2) = max(rindsrc(0:2), 0)
    !     else if (source == 4) then
    !         rloc(0) = xarr(ix) + 0.5_dp * dx
    !         rloc(1) = yarr(iy) + 0.5_dp * dy
    !         rloc(2) = maxz

    !         rindsrc(0) = mod(int( rloc(0) / dx ), nx)
    !         rindsrc(1) = mod(int( rloc(1) / dy ), ny)
    !         rindsrc(2) = mod(int((rloc(2) - 1.0e-8_dp) / dz ), nz)
    !         rind(0:2) = rindsrc(0:2)
    !     end if
    ! end subroutine photon_initloc

    ! subroutine photon_initdir(source, dirsol, dirview, rdir, rdir_sign)
    !     implicit none
    !     integer, intent(in) :: source
    !     real(dp), intent(in) :: dirsol(0:2)
    !     real(dp), intent(in) :: dirview(0:2)
    !     real(dp), intent(out) :: rdir(0:2)
    !     integer, intent(out) :: rdir_sign(0:2)
    !     real(dp) :: rnd, mu, phi, sint, pi
    !     pi = 4.0_dp * atan(1.0_dp)

    !     if (source == 0) then 
    !         rdir = dirsol
    !     else if (source == 1) then
    !         ! Lambertian source
    !         call random_number(rnd)
    !         mu = sqrt(rnd)
    !         call random_number(rnd)
    !         phi = 2.0_dp * pi * rnd
    !         sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
    !         rdir(0) = sint * cos(phi)
    !         rdir(1) = sint * sin(phi)
    !         rdir(2) = -mu
    !     else if (source == 2) then
    !         ! Isotropic source
    !         call random_number(rnd)
    !         mu = 2.0_dp * rnd - 1.0_dp
    !         call random_number(rnd)
    !         phi = 2.0_dp * pi * rnd
    !         sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
    !         rdir(0) = sint * cos(phi)
    !         rdir(1) = sint * sin(phi)
    !         rdir(2) = mu
    !     else if (source == 3) then
    !         ! Lambertian source
    !         call random_number(rnd)
    !         mu = sqrt(rnd)
    !         call random_number(rnd)
    !         phi = 2.0_dp * pi * rnd
    !         sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
    !         if (ia == 0) then ! x+
    !             rdir(0) = mu
    !             rdir(1) = sint * sin(phi)
    !             rdir(2) = sint * cos(phi)
    !         else if (ia == 1) then ! x-
    !             rdir(0) = -mu
    !             rdir(1) = sint * sin(phi)
    !             rdir(2) = -sint * cos(phi)
    !         else if (ia == 2) then ! y+
    !             rdir(0) = sint * cos(phi)
    !             rdir(1) = mu
    !             rdir(2) = sint * sin(phi)
    !         else if (ia == 3) then ! y-
    !             rdir(0) = sint * cos(phi)
    !             rdir(1) = -mu
    !             rdir(2) = -sint * sin(phi)
    !         else if (ia == 4) then ! z+
    !             rdir(0) = sint * cos(phi)
    !             rdir(1) = sint * sin(phi)
    !             rdir(2) = mu
    !         else if (ia == 5) then ! z-
    !             rdir(0) = sint * cos(phi)
    !             rdir(1) = sint * sin(phi)
    !             rdir(2) = -mu
    !         end if
    !     else if (source == 4) then
    !         rdir = dirview
    !     end if
    !     rdir = rdir / sqrt(sum(rdir**2))
    !     rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))
    ! end subroutine photon_initdir

    ! subroutine photon_intersect(nx, ny, nz, r5ind, r5loc, r5dir, r5dir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)
    !     implicit none
    !     integer, intent(in) :: r5ind(0:2)
    !     real(dp), intent(in) :: r5loc(0:2)
    !     real(dp), intent(in) :: r5dir(0:2)
    !     integer, intent(in) :: r5dir_sign(0:2)
    !     real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
    !     real(dp), intent(out) :: rdist
    !     real(dp), intent(out) :: rdloc(0:2)
    !     integer, intent(out) :: icase, isign 
    !     real(dp) :: xbnd(0:1), ybnd(0:1), zbnd(0:1)
    !     real(dp) :: dc(0:2)

    !     xbnd(0) = xarr(r5ind(0))
    !     xbnd(1) = xarr(r5ind(0)+1)
    !     ybnd(0) = yarr(r5ind(1))
    !     ybnd(1) = yarr(r5ind(1)+1)
    !     zbnd(0) = zarr(r5ind(2))
    !     zbnd(1) = zarr(r5ind(2)+1)

    !     dc(0) = xbnd(r5dir_sign(0)) - r5loc(0)
    !     dc(1) = ybnd(r5dir_sign(1)) - r5loc(1)
    !     dc(2) = zbnd(r5dir_sign(2)) - r5loc(2)

    !     icase = 2
    !     icase = merge(1, icase, abs(dc(1) * r5dir(icase)) < abs(dc(icase) * r5dir(1)))
    !     icase = merge(0, icase, abs(dc(0) * r5dir(icase)) < abs(dc(icase) * r5dir(0)))

    !     isign = r5dir_sign(icase)

    !     rdloc(0:2) = dc(icase) * r5dir(0:2) / r5dir(icase)
    !     rdist = sqrt(sum(rdloc**2))
    ! end subroutine photon_intersect
    ! subroutine photon_scattering(g_grid, rdir, rdir_sign, ptau, scaord)
    !     implicit none
    !     real(dp), intent(in) :: g_grid
    !     real(dp), intent(inout) :: rdir(0:2)
    !     integer, intent(inout) :: rdir_sign(0:2)
    !     real(dp), intent(out) :: ptau
    !     integer, intent(inout) :: scaord

    !     real(dp) :: rnd, s, mu, phi, sint, ux, uy, uz, denom
    !     real(dp) :: pi
    !     pi = 4.0_dp * atan(1.0_dp)
    !     ! Henyey-Greenstein
    !     g = max(g_grid, 1.0e-8_dp)
    !     call random_number(rnd); s = 2.0_dp * rnd - 1.0_dp
    !     mu = (1.0_dp + g**2.0_dp - ((1.0_dp - g**2.0_dp) / (1.0_dp + g * s))**2.0_dp) / (2.0_dp * g)
    !     call random_number(rnd); phi = 2.0_dp * pi * rnd
    !     sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
    !     ux = rdir(0); uy = rdir(1); uz = rdir(2)
    !     denom = sqrt(ux*ux + uy*uy)
    !     rdir(0) = merge( &
    !         sint * (-ux * uz * cos(phi) + uy * sin(phi)) / denom + ux * mu, &
    !         merge(sint * cos(phi), -sint * cos(phi), uz > 0.0_dp), &
    !         abs(uz) < 0.9999_dp)
    !     rdir(1) = merge( &
    !         sint * (-uy * uz * cos(phi) - ux * sin(phi)) / denom + uy * mu, &
    !         merge(sint * sin(phi), -sint * sin(phi), uz > 0.0_dp), &
    !         abs(uz) < 0.9999_dp)
    !     rdir(2) = merge( &
    !         sint * cos(phi) * denom + uz * mu, &
    !         merge(mu, -mu, uz > 0.0_dp), &
    !         abs(uz) < 0.9999_dp)
    !     ! normalize
    !     rdir = rdir / sqrt(sum(rdir**2))
    !     rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))
    !     call random_number(rnd); ptau = -log(max(1.0e-300_dp, rnd))
    !     ! weight = new_weight
    !     scaord = scaord + 1
    ! end subroutine photon_scattering

    ! subroutine photon_reflection(rdir, rdir_sign, ptau, scaord)
    !     implicit none
    !     real(dp), intent(inout) :: rdir(0:2)
    !     integer, intent(inout) :: rdir_sign(0:2)
    !     real(dp), intent(out) :: ptau
    !     integer, intent(inout) :: scaord

    !     call random_number(rnd); mu = sqrt(rnd)
    !     call random_number(rnd); phi = 2.0_dp * pi * rnd
    !     sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
    !     rdir(0) = sint * cos(phi)
    !     rdir(1) = sint * sin(phi)
    !     rdir(2) = mu
    !     ! normalize
    !     rdir = rdir / sqrt(sum(rdir**2))
    !     rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))

    !     call random_number(rnd); ptau = -log(max(1.0e-300_dp, rnd))
    !     scaord = scaord + 1
    ! end subroutine photon_reflection

    ! subroutine photon_movegrid(r5ind, r5loc, r5dir, r5dist, icase, isign, nx, ny, dxs, transfer_mode)
    !     implicit none
    !     integer, intent(inout) :: r5ind(0:2)
    !     real(dp), intent(inout) :: r5loc(0:2)
    !     real(dp), intent(in) :: r5dir(0:2)
    !     real(dp), intent(in) :: r5dist
    !     integer, intent(in) :: icase, isign, nx, ny
    !     real(dp), intent(in) :: dxs(0:2)
    !     integer, intent(in) :: transfer_mode
    !     integer :: r5ind2(0:2)
    !     r5loc(0:2) = r5loc(0:2) + r5dir(0:2) * r5dist
    !     r5ind2(0:2) = r5ind(0:2)
    !     r5ind2(icase) = r5ind(icase) + 2 * isign - 1
    !     ! wrap x,y
    !     r5ind(0) = mod(r5ind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * r5ind(0)
    !     r5ind(1) = mod(r5ind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * r5ind(1)
    !     r5ind(2) = r5ind2(2)
    !     ! rloc(0) = mod(rloc(0) + maxx, maxx) * real(isign) + (maxx - mod(maxx - rloc(0), maxx)) * real(1 - isign)
    !     ! rloc(1) = mod(rloc(1) + maxy, maxy) * real(isign) + (maxy - mod(maxy - rloc(1), maxy)) * real(1 - isign)
    !     r5loc(icase) = merge( &
    !         real(r5ind(icase), dp) * dxs(icase) * real(isign, dp) + real(r5ind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
    !         r5loc(icase), &
    !         icase <= 1 )
    ! end subroutine photon_movegrid

    ! subroutine photon_rroulette(new_weight, weight_min, weight_rr, survived, debug)
    !     implicit none
    !     real(dp), intent(inout) :: new_weight
    !     real(dp), intent(in) :: weight_min, weight_rr
    !     logical, intent(out) :: survived
    !     integer, intent(in) :: debug
    !     real(dp) :: rnd
    !     survived = .true.
    !     if (new_weight < weight_min) then
    !         ! Russian Roulette
    !         if (debug == 1) then
    !             print *, "  Russian Roulette triggered. Weight before RR: ", new_weight
    !         end if
    !         call random_number(rnd)
    !         if (rnd * weight_rr > new_weight) then
    !             ! photon terminated
    !             if (debug == 1) then
    !                 print *, "  Photon terminated by RR."
    !             end if
    !             ! exit
    !             survived = .false.
    !         else
    !             new_weight = weight_rr
    !             if (debug == 1) then
    !                 print *, "  Photon survived RR. New Weight: ", new_weight
    !             end if
    !         end if
    !     end if
    ! end subroutine photon_rroulette

    ! subroutine photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
    !     ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
    !     implicit none
    !     real(dp), intent(in) :: rloc(0:2)
    !     integer, intent(in) :: rind(0:2)
    !     real(dp), intent(in) :: dirsol(0:2)
    !     integer, intent(in) :: itlpbmax
    !     real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
    !     real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
    !     integer, intent(in) :: nx, ny, nz
    !     real(dp), intent(in) :: dxs(0:2)
    !     integer, intent(in) :: transfer_mode
    !     integer, intent(in) :: ix, iy, iz, ia, iphoton, it
    !     integer, intent(in) :: debug
    !     integer, intent(in) :: iutraj
    !     real(dp), intent(inout) :: vqllpb
    !     integer :: rlpbind(0:2), rlpbdir_sign(0:2), itlpb, icase, isign
    !     real(dp) :: rlpbloc(0:2), rlpbdir(0:2)
    !     real(dp) :: rdist, rdloc(0:2)
    !     real(dp) :: kext_grid
    !     logical :: survived
    !     vqllpb = 1.0_dp
    !     rlpbloc(0:2) = rloc(0:2)
    !     rlpbdir(0:2) = -dirsol(0:2)
    !     rlpbdir_sign(0:2) = max(0, min(1, ceiling(rlpbdir(0:2))))
    !     rlpbind(0:2) = rind(0:2)
    !     itlpb = 0
    !     do while (itlpb < itlpbmax)
    !         if (debug > 0) then
    !             if (debug == 1) then
    !                 print *, "Local estimation step ", itlpb, " rlpbloc: ", rlpbloc, " rlpbind: ", rlpbind, " rlpbdir: ", rlpbdir
    !             else
    !                 write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'b', itlpb, rlpbloc(0), rlpbloc(1), rlpbloc(2), rlpbind(0), rlpbind(1), rlpbind(2), &
    !                     rlpbdir(0), rlpbdir(1), rlpbdir(2), vqllpb, 0.0_dp
    !             end if
    !         end if
    !         kext_grid = kext(rlpbind(0), rlpbind(1), rlpbind(2))

    !         ! xbnd(0) = xarr(rlpbind(0))
    !         ! xbnd(1) = xarr(rlpbind(0)+1)
    !         ! ybnd(0) = yarr(rlpbind(1))
    !         ! ybnd(1) = yarr(rlpbind(1)+1)
    !         ! zbnd(0) = zarr(rlpbind(2))
    !         ! zbnd(1) = zarr(rlpbind(2)+1)

    !         ! dc(0) = xbnd(rlpbdir_sign(0)) - rlpbloc(0)
    !         ! dc(1) = ybnd(rlpbdir_sign(1)) - rlpbloc(1)
    !         ! dc(2) = zbnd(rlpbdir_sign(2)) - rlpbloc(2)

    !         ! icase = 2
    !         ! icase = merge(1, icase, abs(dc(1) * rlpbdir(icase)) < abs(dc(icase) * rlpbdir(1)))
    !         ! icase = merge(0, icase, abs(dc(0) * rlpbdir(icase)) < abs(dc(icase) * rlpbdir(0)))

    !         ! isign = rlpbdir_sign(icase)

    !         ! rdloc(0:2) = dc(icase) * rlpbdir(0:2) / rlpbdir(icase)
    !         ! rdist = sqrt(sum(rdloc**2))

    !         call photon_intersect(nx, ny, nz, rlpbind, rlpbloc, rlpbdir, rlpbdir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)
    !         ! (r5ind, r5loc, r5dir, r5dir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)


    !         vqllpb = vqllpb * exp(-kext_grid * rdist)

    !         ! ! move to next boundary
    !         ! rlpbloc(0:2) = rlpbloc(0:2) + rlpbdir(0:2) * rdist
    !         ! ! update cell index
    !         ! ! rlpbind(icase) = rlpbind(icase) + isign

    !         ! ! step cell index in icase direction
    !         ! rlpbind2(0:2) = rlpbind(0:2)
    !         ! rlpbind2(icase) = rlpbind(icase) + 2 * isign - 1
    !         ! ! wrap x,y
    !         ! rlpbind(0) = mod(rlpbind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * rlpbind(0)
    !         ! rlpbind(1) = mod(rlpbind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * rlpbind(1)
    !         ! rlpbind(2) = rlpbind2(2)

    !         ! rlpbloc(icase) = merge( &
    !         !     real(rlpbind(icase), dp) * dxs(icase) * real(isign, dp) + real(rlpbind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
    !         !     rlpbloc(icase), &
    !         !     icase <= 1 )
    !         call photon_movegrid(rlpbind, rlpbloc, rlpbdir, rdist, icase, isign, nx, ny, dxs, transfer_mode)

    !         itlpb = itlpb + 1
    !         if (rlpbind(2) >= nz) exit
    !         ! if (vqllpb < 1.0e-4_dp) exit
    !         call photon_rroulette(vqllpb, 0.05_dp, 0.1_dp, survived, debug)
    !         if (.not. survived) exit
    !     end do
    ! end subroutine photon_raytrace

    ! ! subroutine sample_d6(vqla, phflx, rind, idi, iside, nx, ny, nz, swlw)
    ! !     implicit none
    ! !     real(dp), intent(in) :: vqla
    ! !     real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
    ! !     integer, intent(in) :: rind(0:2)
    ! !     integer, intent(in) :: idi, iside
    ! !     integer, intent(in) :: nx, ny, nz, swlw
    ! !     phflx(rind(0), rind(1), rind(2), idi, iside) = phflx(rind(0), rind(1), rind(2), idi, iside) + vqla
    ! ! end subroutine sample_d6

    ! ! subroutine sample_d4(vqla, phconv, rind, nx, ny, nz)
    ! !     implicit none
    ! !     real(dp), intent(in) :: vqla
    ! !     real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
    ! !     integer, intent(in) :: rind(0:2)
    ! !     integer, intent(in) :: nx, ny, nz
    ! !     phconv(rind(0), rind(1), rind(2)) = phconv(rind(0), rind(1), rind(2)) + vqla
    ! ! end subroutine sample_d4

    ! ! subroutine sample_d3(vqla, phimg, rind, nx, ny)
    ! !     implicit none
    ! !     real(dp), intent(in) :: vqla
    ! !     real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
    ! !     integer, intent(in) :: rind(0:2)
    ! !     integer, intent(in) :: nx, ny
    ! !     phimg(rind(0), rind(1)) = phimg(rind(0), rind(1)) + vqla
    ! ! end subroutine sample_d3

    ! ! subroutine sample_scattering(weight_absorbed, new_weight, source, ia, scaord, swlw, &
    ! !     nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
    ! !     rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
    ! !     ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
    ! !     implicit none
    ! !     integer, intent(in) :: source, ia, scaord, swlw
    ! !     integer, intent(in) :: nx, ny, nz
    ! !     integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
    ! !     integer, intent(in) :: transfer_mode
    ! !     integer, intent(in) :: rind(0:2), rindsrc(0:2)
    ! !     real(dp), intent(in) :: weight_absorbed, new_weight
    ! !     real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
    ! !     real(dp), intent(in) :: dxs(0:2)
    ! !     integer, intent(in) :: itlpbmax
    ! !     real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
    ! !     real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
    ! !     real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
    ! !     real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
    ! !     real(dp), intent(inout) :: vqllpb
    ! !     real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
    ! !     real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
    ! !     real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
    ! !     real(dp) :: vqla, notindsrc, mu, g_grid

    ! !     if (source <= 1) then
    ! !         vqla = weight_absorbed
    ! !         call sample_d4(vqla, phconv, rind, nx, ny, nz)
    ! !     else if (source == 2) then
    ! !         notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
    ! !         vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * notindsrc
    ! !         call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)
    ! !     else if (source == 3) then
    ! !         vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * 0.5_dp
    ! !         call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)
    ! !     else if (source == 4) then
    ! !         if (swlw == 0) then
    ! !             vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
    ! !             ! call sample_d3(vqla, phimg, rindsrc, nx, ny)
    ! !         else if (swlw == 1) then
    ! !             call photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
    ! !                 ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
    ! !             mu = -sum(dirsol(0:2) * rdir(0:2))
    ! !             g_grid = gparam(rind(0), rind(1), rind(2))
    ! !             vqla = vqllpb * new_weight * (1.0_dp - g_grid**2) &
    ! !                 / (4.0_dp * pi * (1.0_dp + g_grid**2 - 2.0_dp * g_grid * mu)**1.5_dp * abs(dirsol(2)))
    ! !         end if
    ! !         call sample_d3(vqla, phimg, rindsrc, nx, ny)
    ! !     end if
    ! ! end subroutine sample_scattering

    ! subroutine store_d6(val, recflx, rind, idi, iside, nx, ny, nz, swlw)
    !     implicit none
    !     real(dp), intent(in) :: val(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
    !     real(dp), intent(inout) :: recflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5,0:1)
    !     integer, intent(in) :: rind(0:2)
    !     integer, intent(in) :: idi, iside
    !     integer, intent(in) :: nx, ny, nz, swlw
    !     recflx(rind(0), rind(1), rind(2), idi, iside, 0) = &
    !         recflx(rind(0), rind(1), rind(2), idi, iside, 0) + val(rind(0), rind(1), rind(2), idi, iside)
    !     recflx(rind(0), rind(1), rind(2), idi, iside, 1) = &
    !         recflx(rind(0), rind(1), rind(2), idi, iside, 1) + val(rind(0), rind(1), rind(2), idi, iside)**2.0_dp
    ! end subroutine store_d6

    ! subroutine store_d4(val, recconv, rind, nx, ny, nz)
    !     implicit none
    !     real(dp), intent(in) :: val(0:nx-1,0:ny-1,0:nz-1)
    !     real(dp), intent(inout) :: recconv(0:nx-1,0:ny-1,0:nz-1,0:1)
    !     integer, intent(in) :: rind(0:2)
    !     integer, intent(in) :: nx, ny, nz
    !     recconv(rind(0), rind(1), rind(2), 0) = recconv(rind(0), rind(1), rind(2), 0) + val(rind(0), rind(1), rind(2))
    !     recconv(rind(0), rind(1), rind(2), 1) = recconv(rind(0), rind(1), rind(2), 1) + val(rind(0), rind(1), rind(2))**2.0_dp
    ! end subroutine store_d4

    ! subroutine store_d3(val, recimg, rind, nx, ny)
    !     implicit none
    !     real(dp), intent(in) :: val(0:nx-1,0:ny-1)
    !     real(dp), intent(inout) :: recimg(0:nx-1,0:ny-1,0:1)
    !     integer, intent(in) :: rind(0:2)
    !     integer, intent(in) :: nx, ny
    !     recimg(rind(0), rind(1), 0) = recimg(rind(0), rind(1), 0) + val(rind(0), rind(1))
    !     recimg(rind(0), rind(1), 1) = recimg(rind(0), rind(1), 1) + val(rind(0), rind(1))**2.0_dp
    ! end subroutine store_d3

    ! subroutine record_trajectory(iutraj, ix, iy, iz, ia, iphoton, it, step_type, itlpb, rloc, rind, rdir, weight, ptau)
    !     implicit none
    !     integer, intent(in) :: iutraj
    !     integer, intent(in) :: ix, iy, iz, ia, iphoton, it
    !     character(len=1), intent(in) :: step_type
    !     integer, intent(in) :: itlpb
    !     real(dp), intent(in) :: rloc(0:2)
    !     integer, intent(in) :: rind(0:2)
    !     real(dp), intent(in) :: rdir(0:2)
    !     real(dp), intent(in) :: weight
    !     real(dp), intent(in) :: ptau

    !     write(iutraj,*) ix, iy, iz, ia, iphoton, it, step_type, itlpb, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), &
    !         rdir(0), rdir(1), rdir(2), weight, ptau
    ! end subroutine record_trajectory

    ! subroutine write_output_d6(outrad, nx, ny, nz, nd, na, nphoton, filename)
    !     implicit none
    !     real(dp), intent(in) :: outrad(0:nx-1,0:ny-1,0:nz-1,0:nd-1,0:na-1,0:1)
    !     integer, intent(in) :: nx, ny, nz, nd, na, nphoton
    !     character(len=*), intent(in) :: filename
    !     real(dp) :: val, val2
    !     integer :: unit, ix, iy, iz, idi, ia

    !     unit = 10
    !     open(unit, file = filename, status = "replace", action = "write", form = "formatted")

    !     write(unit, '(A)') "! nx ny nz nd na comp nphoton"
    !     write(unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I10)') nx, ny, nz, nd, na, 2, nphoton
    !     do ix = 0, nx-1
    !         do iy = 0, ny-1
    !             do iz = 0, nz-1
    !                 do idi = 0, nd-1
    !                     do ia = 0, na-1
    !                         val = outrad(ix,iy,iz,idi,ia,0)
    !                         val2 = outrad(ix,iy,iz,idi,ia,1)
    !                         if (abs(val) < 1.0e-20_dp) val = 0.0_dp
    !                         if (abs(val2) < 1.0e-20_dp) val2 = 0.0_dp
    !                         write(unit, '(E15.6, E15.6)') val, val2
    !                     end do
    !                 end do
    !             end do
    !         end do
    !     end do

    !     close(unit)
    ! end subroutine write_output_d6

    ! subroutine write_output_d4(outrad, nx, ny, nz, nphoton, filename)
    !     implicit none
    !     real(dp), intent(in) :: outrad(0:nx-1,0:ny-1,0:nz-1,0:1)
    !     integer, intent(in) :: nx, ny, nz, nphoton
    !     character(len=*), intent(in) :: filename
    !     real(dp) :: val, val2
    !     integer :: unit, ix, iy, iz

    !     unit = 10
    !     open(unit, file = filename, status = "replace", action = "write", form = "formatted")

    !     ! write(unit, '(A)') "! nx ny nz ncase comp nphoton"
    !     ! write(unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I10)') nx, ny, nz, 3, 2, nphoton
    !     write(unit, '(A)') "! nx ny nz comp nphoton"
    !     write(unit, '(I5,1X,I5,1X,I5,1X,I5,1X,I10)') nx, ny, nz, 2, nphoton
    !     do ix = 0, nx-1
    !         do iy = 0, ny-1
    !             do iz = 0, nz-1
    !                 ! do icase = 0, 2
    !                     ! val = outrad(ix,iy,iz,icase,0)
    !                     ! val2 = (outrad(ix,iy,iz,icase,1) - val**2.0_dp) / real(nphoton, dp)
    !                     val = outrad(ix,iy,iz,0)
    !                     val2 = outrad(ix,iy,iz,1)
    !                     if (abs(val) < 1.0e-20_dp) val = 0.0_dp
    !                     if (abs(val2) < 1.0e-20_dp) val2 = 0.0_dp
    !                     write(unit, '(E15.6, E15.6)') val, val2
    !                 ! end do
    !             end do
    !         end do
    !     end do

    !     close(unit)
    ! end subroutine write_output_d4

    ! subroutine write_output_d3(outrad, nx, ny, nphoton, filename)
    !     implicit none
    !     real(dp), intent(in) :: outrad(0:nx-1,0:ny-1,0:1)
    !     integer, intent(in) :: nx, ny, nphoton
    !     character(len=*), intent(in) :: filename
    !     real(dp) :: val, val2
    !     integer :: unit, ix, iy

    !     unit = 10
    !     open(unit, file = filename, status = "replace", action = "write", form = "formatted")

    !     write(unit, '(A)') "! nx ny comp nphoton"
    !     write(unit, '(I5,1X,I5,1X,I5,1X,I10)') nx, ny, 2, nphoton
    !     do ix = 0, nx-1
    !         do iy = 0, ny-1
    !             val = outrad(ix,iy,0)
    !             val2 = outrad(ix,iy,1)
    !             if (abs(val) < 1.0e-20_dp) val = 0.0_dp
    !             if (abs(val2) < 1.0e-20_dp) val2 = 0.0_dp
    !             write(unit, '(E15.6, E15.6)') val, val2
    !         end do
    !     end do

    !     close(unit)
    ! end subroutine write_output_d3

end program test_fullmc




