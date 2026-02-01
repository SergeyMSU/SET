! include "Storage_module.f90"
! include "cuf_kernel.cuf"
! include "cuf_solvers.cuf"
! include "Solvers.f90"

subroutine Print()
    use STORAGE
    implicit none
    
    integer :: i, ierr, unit_num
    real(8) :: x, y, rho, u, v, p
    
    ! Открываем текстовый файл для записи
    open(newunit=unit_num, file='output_data_3.8.txt', status='replace', &
         action='write', iostat=ierr)
    
    if (ierr /= 0) then
        print*, "Ошибка открытия файла для записи"
        return
    endif
    
    ! Записываем заголовок (опционально)
    write(unit_num, '(A)') 'TITLE = "HP"  VARIABLES = x, y, rho, u, v, p, |V|'
    
    ! Проходим по всем ячейкам
    do i = 1, size(host_Cell_par, 2)
        
        ! Извлекаем координаты центра ячейки
        x = host_Cell_center(1, i)  ! x-координата
        y = host_Cell_center(2, i)  ! y-координата
        
        ! Извлекаем физические величины из host_Cell_par
        ! Согласно описанию: индекс 1 - rho, 2 - vx (u), 3 - vy (v), 5 - p
        rho = host_Cell_par(1, i)   ! Плотность
        u   = host_Cell_par(2, i)   ! Компонента скорости по x
        v   = host_Cell_par(3, i)   ! Компонента скорости по y
        p   = host_Cell_par(5, i)   ! Давление
        
        ! Записываем данные в файл
        ! Используем формат с достаточной точностью
        write(unit_num, '(6(ES15.7, 2X))') x, y, rho, u, v, p, sqrt(u**2 + v**2)
        
    end do
    
    ! Закрываем файл
    close(unit_num)
    
    print*, "Данные успешно записаны в файл output_data.txt"
    print*, "Количество записанных ячеек:", size(host_Cell_par, 2)
end subroutine Print


program MIK
    use STORAGE
    use MY_CUDA

    real(8) :: the, dd

    call Set_Storage()
    print*, "Schital"
    ! call Fill_data()

    dd = 1.8

    ! Изменим скорость и плотность на высоких широтах
    if(.TRUE.) then
        do i = 1, size(host_Cell_par(1, :)) 
            if( norm2(host_Cell_center(:, i)) < 0.11) then
                the = polar_angle(host_Cell_center(1, i), host_Cell_center(2, i))

                !if(polar_angle(host_Cell_center(1, i), host_Cell_center(2, i)) > par_pi_8/4.0) then  ! 45
                !if(polar_angle(host_Cell_center(1, i), host_Cell_center(2, i)) > par_pi_8/6.0) then   ! 30
                !if(the > par_pi_8/5.14285714) then   ! 35
                !    host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                !    host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                !    host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2


                ! 27.5 - 32.5
                ! if(the > par_pi_8/5.53846154) then  
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                ! else if(the > par_pi_8/6.5454545454545454) then
                !     dd = 1.0 + 0.8 * (the - par_pi_8/6.5454545454545454) / (par_pi_8 * (1.0/5.53846154 - 1.0/6.5454545454545454))
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                ! end if

                ! 25 - 35
                ! if(the > par_pi_8/5.14285714) then  
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                ! else if(the > par_pi_8/7.2) then
                !     dd = 1.0 + 0.8 * (the - par_pi_8/7.2) / (par_pi_8 * (1.0/5.14285714 - 1.0/7.2))
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                ! end if

                ! 20 - 40
                if(the > par_pi_8/4.5) then  
                    host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                    host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                    host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                else if(the > par_pi_8/9.0) then
                    dd = 1.0 + 0.8 * (the - par_pi_8/7.2) / (par_pi_8 * (1.0/4.5 - 1.0/9.0))
                    host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                    host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                    host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                end if
            end if
        end do
    end if

    call CUDA_info()
    call CUDA_START_MGD()

    call Print()

    call Save_Storage()

end program MIK