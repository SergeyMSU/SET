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
    open(newunit=unit_num, file='output_data_4.1.txt', status='replace', &
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
        x = host_Cell_center(1, i) * 319.319  ! x-координата
        y = host_Cell_center(2, i) * 319.319  ! y-координата
        
        ! Извлекаем физические величины из host_Cell_par
        ! Согласно описанию: индекс 1 - rho, 2 - vx (u), 3 - vy (v), 5 - p
        rho = host_Cell_par(1, i) * 0.073   ! Плотность
        u   = host_Cell_par(2, i) * 10.3565   ! Компонента скорости по x
        v   = host_Cell_par(3, i) * 10.3565   ! Компонента скорости по y
        p   = host_Cell_par(5, i) * 1.30962   ! Давление   10^-13
        
        ! Записываем данные в файл
        ! Используем формат с достаточной точностью
        write(unit_num, '(6(ES15.7, 2X))') x, y, rho, u, v, p, sqrt(u**2 + v**2)
        
    end do
    
    ! Закрываем файл
    close(unit_num)
    
    print*, "Данные успешно записаны в файл output_data.txt"
    print*, "Количество записанных ячеек:", size(host_Cell_par, 2)
end subroutine Print

subroutine Print_mult(step, time)
    use STORAGE
    implicit none
    
    integer, intent(in) :: step
    real(8), intent(in) :: time
    integer :: i, ierr, unit_num, time_unit
    real(8) :: x, y, rho
    character(len=256) :: filename
    character(len=20) :: step_str, time_str

    ! Преобразуем номер шага в строку
    write(step_str, '(I0)') step
    ! Преобразуем время в строку (научный формат)
    write(time_str, '(ES12.5)') time
    time_str = adjustl(time_str)

    ! Имя файла: output_data_3.5_<номер>.txt
    filename = 'output_data_4.1_' // trim(step_str) // '.txt'

    open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
        print*, "Ошибка открытия файла: ", trim(filename)
        return
    endif

    ! Заголовки Tecplot
    write(unit_num, '(A)') 'TITLE = "HP"'
    write(unit_num, '(A)') 'VARIABLES = x, y, rho'
    ! Строка ZONE с временем (и STRANDID=1 на всякий случай, но при одиночной загрузке он не обязателен)
    write(unit_num, '(A, ES12.5)') 'ZONE T="HP", STRANDID=1, SOLUTIONTIME=', time

    print*, "print TIME = ", time
    open(newunit=time_unit, file="times.txt", status='unknown', position='append', action='write', iostat=ierr)
    write(time_unit, '(ES12.5)') time 
    close(time_unit)

    ! Данные
    do i = 1, size(host_Cell_par, 2), 3
        x = host_Cell_center(1, i) * 319.319
        y = host_Cell_center(2, i) * 319.319
        rho = host_Cell_par(1, i) * 0.073
        write(unit_num, '(3(ES15.7, 2X))') x, y, rho
    end do

    close(unit_num)
    ! print*, "Данные записаны в файл ", trim(filename)
    ! print*, "Количество ячеек:", size(host_Cell_par, 2)
end subroutine Print_mult

subroutine write_T_rho(T, rho1, rho2)
    implicit none
    real(8), intent(in) :: T, rho1, rho2          ! Входные аргументы двойной точности
    integer, parameter :: out_unit = 10    ! Номер логического устройства
    character(*), parameter :: filename = "rhoT_4.1.txt"   ! Имя файла

    ! Открыть файл для дозаписи (если не существует – будет создан)
    open(unit=out_unit, file=filename, position='append', action='write', status='unknown')
    
    ! Записать числа в файл в экспоненциальном формате с 15 знаками мантиссы
    write(out_unit, '(3es23.15)') T, rho1, rho2
    
    ! Закрыть файл
    close(out_unit)
end subroutine write_T_rho

program MIK
    use STORAGE
    use MY_CUDA

    real(8) :: the, dd

    call Set_Storage()
    print*, "Schital"
    call Fill_data()

    dd = 1.8

    ! Изменим скорость и плотность на высоких широтах
    if(.FALSE.) then
        do i = 1, size(host_Cell_par(1, :)) 
            if( norm2(host_Cell_center(:, i)) < 0.22) then
                the = polar_angle(host_Cell_center(1, i), host_Cell_center(2, i))

                !if(polar_angle(host_Cell_center(1, i), host_Cell_center(2, i)) > par_pi_8/4.0) then  ! 45
                !if(polar_angle(host_Cell_center(1, i), host_Cell_center(2, i)) > par_pi_8/6.0) then   ! 30
                ! if(the > par_pi_8/5.14285714) then   ! 35
                ! if(the > par_pi_8/4.5) then   ! 40
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                ! end if

                ! 29 - 31
                if(the > par_pi_8/5.80645161) then  
                    host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                    host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                    host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                else if(the > par_pi_8/6.20689655) then
                    dd = 1.0 + 0.8 * (the - par_pi_8/6.20689655) / (par_pi_8 * (1.0/5.80645161 - 1.0/6.20689655))
                    host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                    host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                    host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                end if

                ! ! 27.5 - 32.5
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
                ! if(the > par_pi_8/4.5) then  
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                ! else if(the > par_pi_8/9.0) then
                !     dd = 1.0 + 0.8 * (the - par_pi_8/9.0) / (par_pi_8 * (1.0/4.5 - 1.0/9.0))
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                ! end if

                ! 22.5 - 37.5
                ! if(the > par_pi_8/4.8) then  
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                ! else if(the > par_pi_8/8.0) then
                !     dd = 1.0 + 0.8 * (the - par_pi_8/8.0) / (par_pi_8 * (1.0/4.8 - 1.0/8.0))
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                ! end if

                ! 0 - 30
                ! if(the > par_pi_8/6.0) then  
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * 1.8
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * 1.8
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (1.8)**2
                ! else
                !     dd = 1.0 + 0.8 * (the) / (par_pi_8 /6.0)
                !     host_Cell_par(2, i) = host_Cell_par(2, i) * dd
                !     host_Cell_par(3, i) = host_Cell_par(3, i) * dd
                !     host_Cell_par(1, i) = host_Cell_par(1, i) / (dd)**2
                ! end if

            end if
        end do
    end if

    print*, "CUDA_PROVERKA"
    call flush(6)
    call CUDA_info()
    print*, "CUDA_START_MGD"
    call flush(6)
    call CUDA_START_MGD()

    !call Print()

    !call Save_Storage()

end program MIK