module fullmc_funcs
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: pi = acos(-1.0_dp)

    abstract interface
        subroutine rec_scat_iface(weight_absorbed, new_weight, ia, scaord, swlw, &
            nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
            rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
            ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
            integer, parameter :: dp = selected_real_kind(15, 307)
            integer, intent(in) :: ia, scaord, swlw
            integer, intent(in) :: nx, ny, nz
            integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
            integer, intent(in) :: transfer_mode
            integer, intent(in) :: rind(0:2), rindsrc(0:2)
            real(dp), intent(in) :: weight_absorbed, new_weight
            real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
            real(dp), intent(in) :: dxs(0:2)
            integer, intent(in) :: itlpbmax
            real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
            real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
            real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
            real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
            real(dp), intent(inout) :: vqllpb
            real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
            real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
            real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
        end subroutine rec_scat_iface

        subroutine rec_bound_iface(weight_absorbed, new_weight, ia, scaord, swlw, icase, isign, rdir, &
            rind, rindsrc, nx, ny, nz, phflx, phconv, bplnk)
            implicit none
            integer, parameter :: dp = selected_real_kind(15, 307)
            integer, intent(in) :: ia, scaord, swlw, icase, isign
            integer, intent(in) :: nx, ny, nz
            integer, intent(in) :: rind(0:2), rindsrc(0:2)
            real(dp), intent(in) :: weight_absorbed, new_weight
            real(dp), intent(in) :: rdir(0:2)
            real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
            real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
            real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
            real(dp) :: notindsrc, vqla
        end subroutine rec_bound_iface
    end interface

    contains

    subroutine recorder_scat1(weight_absorbed, new_weight, ia, scaord, swlw, &
            nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
            rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
            ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
        implicit none
        integer, intent(in) :: ia, scaord, swlw
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
        integer, intent(in) :: transfer_mode
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
        real(dp), intent(in) :: dxs(0:2)
        integer, intent(in) :: itlpbmax
        real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        real(dp), intent(inout) :: vqllpb
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
        real(dp) :: notindsrc, vqla, mu, g_grid

        vqla = weight_absorbed
        call sample_d4(vqla, phconv, rind, nx, ny, nz)

    end subroutine recorder_scat1

    subroutine recorder_scat2(weight_absorbed, new_weight, ia, scaord, swlw, &
            nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
            rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
            ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
        implicit none
        integer, intent(in) :: ia, scaord, swlw
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
        integer, intent(in) :: transfer_mode
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
        real(dp), intent(in) :: dxs(0:2)
        integer, intent(in) :: itlpbmax
        real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        real(dp), intent(inout) :: vqllpb
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
        real(dp) :: notindsrc, vqla, mu, g_grid

        notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
        vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * notindsrc
        call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)

    end subroutine recorder_scat2

    subroutine recorder_scat3(weight_absorbed, new_weight, ia, scaord, swlw, &
            nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
            rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
            ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
        implicit none
        integer, intent(in) :: ia, scaord, swlw
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
        integer, intent(in) :: transfer_mode
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
        real(dp), intent(in) :: dxs(0:2)
        integer, intent(in) :: itlpbmax
        real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        real(dp), intent(inout) :: vqllpb
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
        real(dp) :: notindsrc, vqla, mu, g_grid

        vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * 0.5_dp
        call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)

    end subroutine recorder_scat3
    
    subroutine recorder_scat4(weight_absorbed, new_weight, ia, scaord, swlw, &
            nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
            rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
            ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
        implicit none
        integer, intent(in) :: ia, scaord, swlw
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
        integer, intent(in) :: transfer_mode
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
        real(dp), intent(in) :: dxs(0:2)
        integer, intent(in) :: itlpbmax
        real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        real(dp), intent(inout) :: vqllpb
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
        real(dp) :: notindsrc, vqla, mu, g_grid

        if (swlw == 0) then
            vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
            ! call sample_d3(vqla, phimg, rindsrc, nx, ny)
        else if (swlw == 1) then
            call photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
                ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
            mu = -sum(dirsol(0:2) * rdir(0:2))
            g_grid = gparam(rind(0), rind(1), rind(2))
            vqla = vqllpb * new_weight * (1.0_dp - g_grid**2) &
                / (4.0_dp * pi * (1.0_dp + g_grid**2 - 2.0_dp * g_grid * mu)**1.5_dp * abs(dirsol(2)))
        end if
        call sample_d3(vqla, phimg, rindsrc, nx, ny)

    end subroutine recorder_scat4

    subroutine recorder_bound1(weight_absorbed, new_weight, ia, scaord, swlw, icase, isign, rdir, &
            rind, rindsrc, nx, ny, nz, phflx, phconv, bplnk)
        implicit none
        integer, intent(in) :: ia, scaord, swlw, icase, isign
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rdir(0:2)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp) :: notindsrc, vqla
        vqla = new_weight !* abs(rdir(icase))
        call sample_d6(vqla, phflx, rind, min(min(scaord, 1), swlw), 2 * icase + isign, nx, ny, nz, swlw)
        vqla = weight_absorbed
        call sample_d4(vqla, phconv, rind, nx, ny, nz)
    end subroutine recorder_bound1

    subroutine recorder_bound2(weight_absorbed, new_weight, ia, scaord, swlw, icase, isign, rdir, &
            rind, rindsrc, nx, ny, nz, phflx, phconv, bplnk)
        implicit none
        integer, intent(in) :: ia, scaord, swlw, icase, isign
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rdir(0:2)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp) :: notindsrc, vqla
        notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
        vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * notindsrc
        call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)
    end subroutine recorder_bound2

    subroutine recorder_bound3(weight_absorbed, new_weight, ia, scaord, swlw, icase, isign, rdir, &
            rind, rindsrc, nx, ny, nz, phflx, phconv, bplnk)
        implicit none
        integer, intent(in) :: ia, scaord, swlw, icase, isign
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rdir(0:2)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp) :: notindsrc, vqla
        vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * 0.5_dp
        call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)
    end subroutine recorder_bound3

    subroutine recorder_bound4(weight_absorbed, new_weight, ia, scaord, swlw, icase, isign, rdir, &
            rind, rindsrc, nx, ny, nz, phflx, phconv, bplnk)
        implicit none
        integer, intent(in) :: ia, scaord, swlw, icase, isign
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: rind(0:2), rindsrc(0:2)
        real(dp), intent(in) :: weight_absorbed, new_weight
        real(dp), intent(in) :: rdir(0:2)
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
        real(dp) :: notindsrc, vqla
        ! Do nothing...
    end subroutine recorder_bound4

    subroutine sample_d6(vqla, phflx, rind, idi, iside, nx, ny, nz, swlw)
        implicit none
        real(dp), intent(in) :: vqla
        real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        integer, intent(in) :: rind(0:2)
        integer, intent(in) :: idi, iside
        integer, intent(in) :: nx, ny, nz, swlw
        phflx(rind(0), rind(1), rind(2), idi, iside) = phflx(rind(0), rind(1), rind(2), idi, iside) + vqla
    end subroutine sample_d6

    subroutine sample_d4(vqla, phconv, rind, nx, ny, nz)
        implicit none
        real(dp), intent(in) :: vqla
        real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
        integer, intent(in) :: rind(0:2)
        integer, intent(in) :: nx, ny, nz
        phconv(rind(0), rind(1), rind(2)) = phconv(rind(0), rind(1), rind(2)) + vqla
    end subroutine sample_d4

    subroutine sample_d3(vqla, phimg, rind, nx, ny)
        implicit none
        real(dp), intent(in) :: vqla
        real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
        integer, intent(in) :: rind(0:2)
        integer, intent(in) :: nx, ny
        phimg(rind(0), rind(1)) = phimg(rind(0), rind(1)) + vqla
    end subroutine sample_d3

        subroutine photon_initloc(source, ix, iy, iz, ia, xarr, yarr, zarr, dx, dy, dz, nx, ny, nz, maxx, maxy, maxz, rloc, rindsrc, rind)
        implicit none
        integer, intent(in) :: source, ix, iy, iz, ia, nx, ny, nz
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        real(dp), intent(in) :: dx, dy, dz
        real(dp), intent(in) :: maxx, maxy, maxz
        real(dp), intent(out) :: rloc(0:2)
        integer, intent(out) :: rindsrc(0:2), rind(0:2)
        real(dp) :: rnd
        if (source <= 1) then
            ! Top of model
            call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
            call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
            rloc(0) = mod(rloc(0) + maxx, maxx)
            rloc(1) = mod(rloc(1) + maxy, maxy)
            rloc(2) = maxz

            rindsrc(0) = mod(int( rloc(0) / dx ), nx)
            rindsrc(1) = mod(int( rloc(1) / dy ), ny)
            rindsrc(2) = nz
            rind(0) = rindsrc(0)
            rind(1) = rindsrc(1)
            rind(2) = rindsrc(2) - 1
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
        else if (source == 3) then
            ! Grid sides
            if (ia == 0) then ! x+
                rloc(0) = xarr(ix)
                call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
                call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
                rindsrc(0) = mod(int((rloc(0) + 1.0e-8_dp) / dx ), nx)
                rindsrc(1) = mod(int( rloc(1) / dy ), ny)
                rindsrc(2) = mod(int( rloc(2) / dz ), nz)
            else if (ia == 1) then ! x-
                rloc(0) = xarr(ix + 1)
                call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
                call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
                rindsrc(0) = mod(int((rloc(0) - 1.0e-8_dp) / dx ), nx)
                rindsrc(1) = mod(int( rloc(1) / dy ), ny)
                rindsrc(2) = mod(int( rloc(2) / dz ), nz)
            else if (ia == 2) then ! y+
                call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
                rloc(1) = yarr(iy)
                call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
                rindsrc(0) = mod(int( rloc(0) / dx ), nx)
                rindsrc(1) = mod(int((rloc(1) + 1.0e-8_dp) / dy ), ny)
                rindsrc(2) = mod(int( rloc(2) / dz ), nz)
            else if (ia == 3) then ! y-
                call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
                rloc(1) = yarr(iy + 1)
                call random_number(rnd); rloc(2) = zarr(iz) + rnd * dz
                rindsrc(0) = mod(int( rloc(0) / dx ), nx)
                rindsrc(1) = mod(int((rloc(1) - 1.0e-8_dp) / dy ), ny)
                rindsrc(2) = mod(int( rloc(2) / dz ), nz)
            else if (ia == 4) then ! z+
                call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
                call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
                rloc(2) = zarr(iz)
                rindsrc(0) = mod(int( rloc(0) / dx ), nx)
                rindsrc(1) = mod(int( rloc(1) / dy ), ny)
                rindsrc(2) = int(floor((rloc(2) + 1.0e-8_dp) / dz ))
            else if (ia == 5) then ! z-
                call random_number(rnd); rloc(0) = xarr(ix) + rnd * dx
                call random_number(rnd); rloc(1) = yarr(iy) + rnd * dy
                rloc(2) = zarr(iz + 1)
                rindsrc(0) = mod(int( rloc(0) / dx ), nx)
                rindsrc(1) = mod(int( rloc(1) / dy ), ny)
                rindsrc(2) = int(floor((rloc(2) - 1.0e-8_dp) / dz ))
            end if
            rind(0:2) = max(rindsrc(0:2), 0)
        else if (source == 4) then
            rloc(0) = xarr(ix) + 0.5_dp * dx
            rloc(1) = yarr(iy) + 0.5_dp * dy
            rloc(2) = maxz

            rindsrc(0) = mod(int( rloc(0) / dx ), nx)
            rindsrc(1) = mod(int( rloc(1) / dy ), ny)
            rindsrc(2) = mod(int((rloc(2) - 1.0e-8_dp) / dz ), nz)
            rind(0:2) = rindsrc(0:2)
        end if
    end subroutine photon_initloc

    subroutine photon_initdir(source, dirsol, dirview, rdir, rdir_sign, ia)
        implicit none
        integer, intent(in) :: source
        real(dp), intent(in) :: dirsol(0:2)
        real(dp), intent(in) :: dirview(0:2)
        integer, intent(in) :: ia
        real(dp), intent(out) :: rdir(0:2)
        integer, intent(out) :: rdir_sign(0:2)
        real(dp) :: rnd, mu, phi, sint, pi
        pi = 4.0_dp * atan(1.0_dp)

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
        else if (source == 3) then
            ! Lambertian source
            call random_number(rnd)
            mu = sqrt(rnd)
            call random_number(rnd)
            phi = 2.0_dp * pi * rnd
            sint = sqrt(max(0.0_dp, 1.0_dp - mu**2))
            if (ia == 0) then ! x+
                rdir(0) = mu
                rdir(1) = sint * sin(phi)
                rdir(2) = sint * cos(phi)
            else if (ia == 1) then ! x-
                rdir(0) = -mu
                rdir(1) = sint * sin(phi)
                rdir(2) = -sint * cos(phi)
            else if (ia == 2) then ! y+
                rdir(0) = sint * cos(phi)
                rdir(1) = mu
                rdir(2) = sint * sin(phi)
            else if (ia == 3) then ! y-
                rdir(0) = sint * cos(phi)
                rdir(1) = -mu
                rdir(2) = -sint * sin(phi)
            else if (ia == 4) then ! z+
                rdir(0) = sint * cos(phi)
                rdir(1) = sint * sin(phi)
                rdir(2) = mu
            else if (ia == 5) then ! z-
                rdir(0) = sint * cos(phi)
                rdir(1) = sint * sin(phi)
                rdir(2) = -mu
            end if
        else if (source == 4) then
            rdir = dirview
        end if
        rdir = rdir / sqrt(sum(rdir**2))
        rdir_sign(0:2) = max(0, min(1, ceiling(rdir(0:2))))
    end subroutine photon_initdir

    subroutine photon_intersect(nx, ny, nz, r5ind, r5loc, r5dir, r5dir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)
        implicit none
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: r5ind(0:2)
        real(dp), intent(in) :: r5loc(0:2)
        real(dp), intent(in) :: r5dir(0:2)
        integer, intent(in) :: r5dir_sign(0:2)
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        real(dp), intent(out) :: rdist
        real(dp), intent(out) :: rdloc(0:2)
        integer, intent(out) :: icase, isign 
        real(dp) :: xbnd(0:1), ybnd(0:1), zbnd(0:1)
        real(dp) :: dc(0:2)

        xbnd(0) = xarr(r5ind(0))
        xbnd(1) = xarr(r5ind(0)+1)
        ybnd(0) = yarr(r5ind(1))
        ybnd(1) = yarr(r5ind(1)+1)
        zbnd(0) = zarr(r5ind(2))
        zbnd(1) = zarr(r5ind(2)+1)

        dc(0) = xbnd(r5dir_sign(0)) - r5loc(0)
        dc(1) = ybnd(r5dir_sign(1)) - r5loc(1)
        dc(2) = zbnd(r5dir_sign(2)) - r5loc(2)

        icase = 2
        icase = merge(1, icase, abs(dc(1) * r5dir(icase)) < abs(dc(icase) * r5dir(1)))
        icase = merge(0, icase, abs(dc(0) * r5dir(icase)) < abs(dc(icase) * r5dir(0)))

        isign = r5dir_sign(icase)

        rdloc(0:2) = dc(icase) * r5dir(0:2) / r5dir(icase)
        rdist = sqrt(sum(rdloc**2))
    end subroutine photon_intersect

    subroutine photon_scattering(g_grid, rdir, rdir_sign, ptau, scaord)
        implicit none
        real(dp), intent(in) :: g_grid
        real(dp), intent(inout) :: rdir(0:2)
        integer, intent(inout) :: rdir_sign(0:2)
        real(dp), intent(out) :: ptau
        integer, intent(inout) :: scaord

        real(dp) :: rnd, s, mu, phi, sint, ux, uy, uz, denom
        real(dp) :: g

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
        ! weight = new_weight
        scaord = scaord + 1
    end subroutine photon_scattering

    subroutine photon_reflection(rdir, rdir_sign, ptau, scaord)
        implicit none
        real(dp), intent(inout) :: rdir(0:2)
        integer, intent(inout) :: rdir_sign(0:2)
        real(dp), intent(out) :: ptau
        integer, intent(inout) :: scaord
        real(dp) :: rnd, mu, phi, sint

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
    end subroutine photon_reflection

    subroutine photon_movegrid(r5ind, r5loc, r5dir, r5dist, icase, isign, nx, ny, dxs, transfer_mode)
        implicit none
        integer, intent(inout) :: r5ind(0:2)
        real(dp), intent(inout) :: r5loc(0:2)
        real(dp), intent(in) :: r5dir(0:2)
        real(dp), intent(in) :: r5dist
        integer, intent(in) :: icase, isign, nx, ny
        real(dp), intent(in) :: dxs(0:2)
        integer, intent(in) :: transfer_mode
        integer :: r5ind2(0:2)
        r5loc(0:2) = r5loc(0:2) + r5dir(0:2) * r5dist
        r5ind2(0:2) = r5ind(0:2)
        r5ind2(icase) = r5ind(icase) + 2 * isign - 1
        ! wrap x,y
        r5ind(0) = mod(r5ind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * r5ind(0)
        r5ind(1) = mod(r5ind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * r5ind(1)
        r5ind(2) = r5ind2(2)
        ! rloc(0) = mod(rloc(0) + maxx, maxx) * real(isign) + (maxx - mod(maxx - rloc(0), maxx)) * real(1 - isign)
        ! rloc(1) = mod(rloc(1) + maxy, maxy) * real(isign) + (maxy - mod(maxy - rloc(1), maxy)) * real(1 - isign)
        r5loc(icase) = merge( &
            real(r5ind(icase), dp) * dxs(icase) * real(isign, dp) + real(r5ind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
            r5loc(icase), &
            icase <= 1 )
    end subroutine photon_movegrid

    subroutine photon_rroulette(new_weight, weight_min, weight_rr, survived, debug)
        implicit none
        real(dp), intent(inout) :: new_weight
        real(dp), intent(in) :: weight_min, weight_rr
        logical, intent(out) :: survived
        integer, intent(in) :: debug
        real(dp) :: rnd
        survived = .true.
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
                ! exit
                survived = .false.
            else
                new_weight = weight_rr
                if (debug == 1) then
                    print *, "  Photon survived RR. New Weight: ", new_weight
                end if
            end if
        end if
    end subroutine photon_rroulette

    subroutine photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
        ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
        implicit none
        real(dp), intent(in) :: rloc(0:2)
        integer, intent(in) :: rind(0:2)
        real(dp), intent(in) :: dirsol(0:2)
        integer, intent(in) :: itlpbmax
        real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(in) :: dxs(0:2)
        integer, intent(in) :: transfer_mode
        integer, intent(in) :: ix, iy, iz, ia, iphoton, it
        integer, intent(in) :: debug
        integer, intent(in) :: iutraj
        real(dp), intent(inout) :: vqllpb
        integer :: rlpbind(0:2), rlpbdir_sign(0:2), itlpb, icase, isign
        real(dp) :: rlpbloc(0:2), rlpbdir(0:2)
        real(dp) :: rdist, rdloc(0:2)
        real(dp) :: kext_grid
        logical :: survived
        vqllpb = 1.0_dp
        rlpbloc(0:2) = rloc(0:2)
        rlpbdir(0:2) = -dirsol(0:2)
        rlpbdir_sign(0:2) = max(0, min(1, ceiling(rlpbdir(0:2))))
        rlpbind(0:2) = rind(0:2)
        itlpb = 0
        do while (itlpb < itlpbmax)
            if (debug > 0) then
                if (debug == 1) then
                    print *, "Local estimation step ", itlpb, " rlpbloc: ", rlpbloc, " rlpbind: ", rlpbind, " rlpbdir: ", rlpbdir
                else
                    write(iutraj,*) ix, iy, iz, ia, iphoton, it, 'b', itlpb, rlpbloc(0), rlpbloc(1), rlpbloc(2), rlpbind(0), rlpbind(1), rlpbind(2), &
                        rlpbdir(0), rlpbdir(1), rlpbdir(2), vqllpb, 0.0_dp
                end if
            end if
            kext_grid = kext(rlpbind(0), rlpbind(1), rlpbind(2))

            ! xbnd(0) = xarr(rlpbind(0))
            ! xbnd(1) = xarr(rlpbind(0)+1)
            ! ybnd(0) = yarr(rlpbind(1))
            ! ybnd(1) = yarr(rlpbind(1)+1)
            ! zbnd(0) = zarr(rlpbind(2))
            ! zbnd(1) = zarr(rlpbind(2)+1)

            ! dc(0) = xbnd(rlpbdir_sign(0)) - rlpbloc(0)
            ! dc(1) = ybnd(rlpbdir_sign(1)) - rlpbloc(1)
            ! dc(2) = zbnd(rlpbdir_sign(2)) - rlpbloc(2)

            ! icase = 2
            ! icase = merge(1, icase, abs(dc(1) * rlpbdir(icase)) < abs(dc(icase) * rlpbdir(1)))
            ! icase = merge(0, icase, abs(dc(0) * rlpbdir(icase)) < abs(dc(icase) * rlpbdir(0)))

            ! isign = rlpbdir_sign(icase)

            ! rdloc(0:2) = dc(icase) * rlpbdir(0:2) / rlpbdir(icase)
            ! rdist = sqrt(sum(rdloc**2))

            call photon_intersect(nx, ny, nz, rlpbind, rlpbloc, rlpbdir, rlpbdir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)
            ! (r5ind, r5loc, r5dir, r5dir_sign, xarr, yarr, zarr, rdist, rdloc, icase, isign)


            vqllpb = vqllpb * exp(-kext_grid * rdist)

            ! ! move to next boundary
            ! rlpbloc(0:2) = rlpbloc(0:2) + rlpbdir(0:2) * rdist
            ! ! update cell index
            ! ! rlpbind(icase) = rlpbind(icase) + isign

            ! ! step cell index in icase direction
            ! rlpbind2(0:2) = rlpbind(0:2)
            ! rlpbind2(icase) = rlpbind(icase) + 2 * isign - 1
            ! ! wrap x,y
            ! rlpbind(0) = mod(rlpbind2(0) + nx, nx) * transfer_mode + (1 - transfer_mode) * rlpbind(0)
            ! rlpbind(1) = mod(rlpbind2(1) + ny, ny) * transfer_mode + (1 - transfer_mode) * rlpbind(1)
            ! rlpbind(2) = rlpbind2(2)

            ! rlpbloc(icase) = merge( &
            !     real(rlpbind(icase), dp) * dxs(icase) * real(isign, dp) + real(rlpbind(icase) + 1, dp) * dxs(icase) * real(1 - isign, dp), &
            !     rlpbloc(icase), &
            !     icase <= 1 )
            call photon_movegrid(rlpbind, rlpbloc, rlpbdir, rdist, icase, isign, nx, ny, dxs, transfer_mode)

            itlpb = itlpb + 1
            if (rlpbind(2) >= nz) exit
            ! if (vqllpb < 1.0e-4_dp) exit
            call photon_rroulette(vqllpb, 0.05_dp, 0.1_dp, survived, debug)
            if (.not. survived) exit
        end do
    end subroutine photon_raytrace

    ! subroutine sample_d6(vqla, phflx, rind, idi, iside, nx, ny, nz, swlw)
    !     implicit none
    !     real(dp), intent(in) :: vqla
    !     real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
    !     integer, intent(in) :: rind(0:2)
    !     integer, intent(in) :: idi, iside
    !     integer, intent(in) :: nx, ny, nz, swlw
    !     phflx(rind(0), rind(1), rind(2), idi, iside) = phflx(rind(0), rind(1), rind(2), idi, iside) + vqla
    ! end subroutine sample_d6

    ! subroutine sample_d4(vqla, phconv, rind, nx, ny, nz)
    !     implicit none
    !     real(dp), intent(in) :: vqla
    !     real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
    !     integer, intent(in) :: rind(0:2)
    !     integer, intent(in) :: nx, ny, nz
    !     phconv(rind(0), rind(1), rind(2)) = phconv(rind(0), rind(1), rind(2)) + vqla
    ! end subroutine sample_d4

    ! subroutine sample_d3(vqla, phimg, rind, nx, ny)
    !     implicit none
    !     real(dp), intent(in) :: vqla
    !     real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
    !     integer, intent(in) :: rind(0:2)
    !     integer, intent(in) :: nx, ny
    !     phimg(rind(0), rind(1)) = phimg(rind(0), rind(1)) + vqla
    ! end subroutine sample_d3

    ! subroutine sample_scattering(weight_absorbed, new_weight, source, ia, scaord, swlw, &
    !     nx, ny, nz, rind, rindsrc, kext, bplnk, gparam, phconv, phflx, phimg, &
    !     rloc, rdir, dirsol, itlpbmax, xarr, yarr, zarr, dxs, transfer_mode, &
    !     ix, iy, iz, iphoton, it, debug, iutraj, vqllpb)
    !     implicit none
    !     integer, intent(in) :: source, ia, scaord, swlw
    !     integer, intent(in) :: nx, ny, nz
    !     integer, intent(in) :: ix, iy, iz, iphoton, it, debug, iutraj
    !     integer, intent(in) :: transfer_mode
    !     integer, intent(in) :: rind(0:2), rindsrc(0:2)
    !     real(dp), intent(in) :: weight_absorbed, new_weight
    !     real(dp), intent(in) :: rloc(0:2), rdir(0:2), dirsol(0:2)
    !     real(dp), intent(in) :: dxs(0:2)
    !     integer, intent(in) :: itlpbmax
    !     real(dp), intent(in) :: kext(0:nx-1,0:ny-1,0:nz-1)
    !     real(dp), intent(in) :: bplnk(0:nx-1,0:ny-1,0:nz-1)
    !     real(dp), intent(in) :: gparam(0:nx-1,0:ny-1,0:nz-1)
    !     real(dp), intent(in) :: xarr(0:nx), yarr(0:ny), zarr(0:nz)
    !     real(dp), intent(inout) :: vqllpb
    !     real(dp), intent(inout) :: phconv(0:nx-1,0:ny-1,0:nz-1)
    !     real(dp), intent(inout) :: phflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
    !     real(dp), intent(inout) :: phimg(0:nx-1,0:ny-1)
    !     real(dp) :: vqla, notindsrc, mu, g_grid

    !     if (source <= 1) then
    !         vqla = weight_absorbed
    !         call sample_d4(vqla, phconv, rind, nx, ny, nz)
    !     else if (source == 2) then
    !         notindsrc = 1.0_dp ! real(min(1, sum(abs(rind(0:2) - rindsrc(0:2)))), dp)
    !         vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * notindsrc
    !         call sample_d4(vqla, phconv, rindsrc, nx, ny, nz)
    !     else if (source == 3) then
    !         vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2)) * 0.5_dp
    !         call sample_d6(vqla, phflx, rindsrc, min(min(scaord, 1), swlw), ia, nx, ny, nz, swlw)
    !     else if (source == 4) then
    !         if (swlw == 0) then
    !             vqla = weight_absorbed * bplnk(rind(0), rind(1), rind(2))
    !             ! call sample_d3(vqla, phimg, rindsrc, nx, ny)
    !         else if (swlw == 1) then
    !             call photon_raytrace(rloc, rind, dirsol, itlpbmax, kext, xarr, yarr, zarr, nx, ny, nz, dxs, transfer_mode, &
    !                 ix, iy, iz, ia, iphoton, it, debug, iutraj, vqllpb)
    !             mu = -sum(dirsol(0:2) * rdir(0:2))
    !             g_grid = gparam(rind(0), rind(1), rind(2))
    !             vqla = vqllpb * new_weight * (1.0_dp - g_grid**2) &
    !                 / (4.0_dp * pi * (1.0_dp + g_grid**2 - 2.0_dp * g_grid * mu)**1.5_dp * abs(dirsol(2)))
    !         end if
    !         call sample_d3(vqla, phimg, rindsrc, nx, ny)
    !     end if
    ! end subroutine sample_scattering

    subroutine store_d6(val, recflx, rind, idi, iside, nx, ny, nz, swlw)
        implicit none
        real(dp), intent(in) :: val(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5)
        real(dp), intent(inout) :: recflx(0:nx-1,0:ny-1,-1:nz,0:swlw,0:5,0:1)
        integer, intent(in) :: rind(0:2)
        integer, intent(in) :: idi, iside
        integer, intent(in) :: nx, ny, nz, swlw
        recflx(rind(0), rind(1), rind(2), idi, iside, 0) = &
            recflx(rind(0), rind(1), rind(2), idi, iside, 0) + val(rind(0), rind(1), rind(2), idi, iside)
        recflx(rind(0), rind(1), rind(2), idi, iside, 1) = &
            recflx(rind(0), rind(1), rind(2), idi, iside, 1) + val(rind(0), rind(1), rind(2), idi, iside)**2.0_dp
    end subroutine store_d6

    subroutine store_d4(val, recconv, rind, nx, ny, nz)
        implicit none
        real(dp), intent(in) :: val(0:nx-1,0:ny-1,0:nz-1)
        real(dp), intent(inout) :: recconv(0:nx-1,0:ny-1,0:nz-1,0:1)
        integer, intent(in) :: rind(0:2)
        integer, intent(in) :: nx, ny, nz
        recconv(rind(0), rind(1), rind(2), 0) = recconv(rind(0), rind(1), rind(2), 0) + val(rind(0), rind(1), rind(2))
        recconv(rind(0), rind(1), rind(2), 1) = recconv(rind(0), rind(1), rind(2), 1) + val(rind(0), rind(1), rind(2))**2.0_dp
    end subroutine store_d4

    subroutine store_d3(val, recimg, rind, nx, ny)
        implicit none
        real(dp), intent(in) :: val(0:nx-1,0:ny-1)
        real(dp), intent(inout) :: recimg(0:nx-1,0:ny-1,0:1)
        integer, intent(in) :: rind(0:2)
        integer, intent(in) :: nx, ny
        recimg(rind(0), rind(1), 0) = recimg(rind(0), rind(1), 0) + val(rind(0), rind(1))
        recimg(rind(0), rind(1), 1) = recimg(rind(0), rind(1), 1) + val(rind(0), rind(1))**2.0_dp
    end subroutine store_d3

    subroutine record_trajectory(iutraj, ix, iy, iz, ia, iphoton, it, step_type, itlpb, rloc, rind, rdir, weight, ptau)
        implicit none
        integer, intent(in) :: iutraj
        integer, intent(in) :: ix, iy, iz, ia, iphoton, it
        character(len=1), intent(in) :: step_type
        integer, intent(in) :: itlpb
        real(dp), intent(in) :: rloc(0:2)
        integer, intent(in) :: rind(0:2)
        real(dp), intent(in) :: rdir(0:2)
        real(dp), intent(in) :: weight
        real(dp), intent(in) :: ptau

        write(iutraj,*) ix, iy, iz, ia, iphoton, it, step_type, itlpb, rloc(0), rloc(1), rloc(2), rind(0), rind(1), rind(2), &
            rdir(0), rdir(1), rdir(2), weight, ptau
    end subroutine record_trajectory

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
                            val2 = outrad(ix,iy,iz,idi,ia,1)
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
                        val2 = outrad(ix,iy,iz,1)
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
end module fullmc_funcs