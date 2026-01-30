! include "Storage_module.f90"
! include "cuf_kernel.cuf"
! include "cuf_solvers.cuf"
! include "Solvers.f90"



program MIK
    use STORAGE
    use MY_CUDA

    call Set_Storage()
    print*, "Schital"
    call Fill_data()

    ! Изменим скорость и плотность на высоких широтах
    if(.FALSE.) then
        do i = 1, size(host_Cell_par(1, :)) 
            if( norm2(host_Cell_center(:, i)) < 0.1) then
                if(polar_angle(host_Cell_center(1, i), host_Cell_center(2, i)) > par_pi_8/4.0) then
                    host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                    host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                    host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                end if
            end if
        end do
    end if

    call CUDA_info()
    call CUDA_START_MGD()

    call Save_Storage()

end program MIK