

module STORAGE 

    real(8), parameter :: par_pi_8 = acos(-1.0_8)         
    real(8), parameter :: par_pi_4 = acos(-1.0_4)      
	real(8), parameter :: par_sqrtpi = sqrt(par_pi_8)
    real(8), parameter :: cpi4 = 12.56637061435917295384_8
    real(8), parameter :: ggg = (5.0_8/3.0_8)
    real(8), parameter :: par_R_character = 1.0_8


    integer(4), parameter  :: host_num_param = 16     ! Число параметров в ячейке


    real(8) :: host_time_all                       ! Текущее время расчёта
    INTEGER(4) :: host_N_cell                       ! Число ячеек в сетке
    INTEGER(4) :: host_N_gran                       ! Число граней в сетке

	real(8), allocatable :: host_Cell_par(:, :)     ! (8, :)

	real(8), allocatable :: host_Cell_center(:, :)      ! (3, :)
	real(8), allocatable :: host_Cell_Volume(:)         ! (:)
    integer(4), allocatable :: host_Cell_gran(:,:)            ! (6, :) Набор из 6 граней для каждой ячейки (если номер = 0, то грани нет в этом направлении)
	real(8), allocatable :: host_Cell_dist(:)                 ! Минимальное расстояние до грани в каждой ячейки 


	real(8), allocatable :: host_Gran_normal(:,:)       ! (3, :) Нормаль грани   
    real(8), allocatable :: host_Gran_square(:)         ! (:) Площадь грани
	real(8), allocatable :: host_Gran_center(:,:)			! (3, :)

	integer(4), allocatable :: host_Gran_neighbour(:,:) ! Соседи каждой грани (2,:) имеют по 2 соседа, нормаль ведёт от первого ко второму
	integer(4), allocatable :: host_Gran_neighbour_TVD(:,:) ! TVD-Соседи каждой грани (2,:) имеют по 2 соседа

	integer(4), allocatable :: host_Gran_type(:)                 ! Показывает тип грани

	real(8), allocatable :: host_Gran_POTOK(:, :)       ! (9, :) поток грани    последний - дивергенция магнитного поля для очистки



    contains

    subroutine Set_Storage()
        real(8) :: vv
        integer(4) :: i
        open(1, file = "FCMHD_1.bin", FORM = 'BINARY', ACTION = "READ")

        read(1) host_time_all
        read(1) host_N_cell
        read(1) host_N_gran

        print*, "AA = ", host_time_all, host_N_cell, host_N_gran

        allocate(host_Cell_par(host_num_param, host_N_cell))
        allocate(host_Cell_center(3, host_N_cell))
        allocate(host_Cell_Volume(host_N_cell))
        allocate(host_Cell_gran(6, host_N_cell))


        allocate(host_Gran_normal(3, host_N_gran))
        allocate(host_Gran_square(host_N_gran))
        allocate(host_Gran_center(3, host_N_gran))
        allocate(host_Gran_neighbour(2, host_N_gran))
        allocate(host_Gran_neighbour_TVD(2, host_N_gran))
        allocate(host_Gran_type(host_N_gran))
        allocate(host_Gran_POTOK(9, host_N_gran))

        read(1) host_Cell_par
        !print*, "S1 = ", host_Cell_par(1, 1), host_Cell_par(2, 1), host_Cell_par(3, 1), host_Cell_par(4, 1), host_Cell_par(5, 1), host_Cell_par(6, 1), host_Cell_par(7, 1), host_Cell_par(8, 1)
        read(1) host_Cell_center
        read(1) host_Cell_Volume
        read(1) host_Cell_gran

        read(1) host_Gran_normal
        read(1) host_Gran_square
        read(1) host_Gran_center
        read(1) host_Gran_neighbour
        read(1) host_Gran_neighbour_TVD
        !print*, "S2 = ", host_Gran_neighbour_TVD(1, 2), host_Gran_neighbour_TVD(2, 2)
        read(1) vv
        !print*, "VV = ", vv
        read(1) host_Gran_type
        read(1) vv
        !print*, "VV = ", vv
        read(1) host_Gran_POTOK

        read(1) vv
        print*, "VV (123)  =  ", vv

        close(1)

        do i = 1, size(host_Cell_par(1, :)) 
            if(host_Cell_par(9, i) < 0.0) then
                print*, "ERROR host_Cell_par 9 < 0  ", host_Cell_par(9, i)
                STOP
            end if
        end do

        do i = 1, size(host_Cell_par(1, :)) 
            if(host_Cell_par(12, i) < 0.0) then
                print*, "ERROR host_Cell_par 12 < 0  ", host_Cell_par(12, i)
                STOP
            end if
        end do

        do i = 1, size(host_Cell_par(1, :)) 
            if(host_Cell_par(16, i) < 0.0) then
                print*, "ERROR host_Cell_par 16 < 0  ", host_Cell_par(16, i)
                STOP
            end if
        end do

        do i = 1, size(host_Cell_par(1, :)) 
            if(host_Cell_par(13, i) < 0.0) then
                print*, "ERROR host_Cell_par 13 < 0  ", host_Cell_par(13, i)
                STOP
            end if
        end do

        do i = 1, size(host_Cell_par(1, :)) 
            if(host_Cell_par(1, i) < 0.000000001) then
                print*, "rho < 0  =  ", host_Cell_par(1, i)
                STOP
            end if
        end do

        do i = 1, size(host_Gran_neighbour(1, :)) 
            if(host_Gran_neighbour(1, i) > host_N_cell) then
                print*, "ERROR host_Gran_neighbour", host_Gran_neighbour(1, i), host_N_cell
                STOP
            end if
        end do

        do i = 1, size(host_Gran_neighbour(1, :)) 
            if(host_Gran_neighbour(1, i) < 1) then
                print*, "ERROR host_Gran_neighbour  ", host_Gran_neighbour(1, i)
                STOP
            end if
        end do

        do i = 1, size(host_Gran_neighbour_TVD(1, :)) 
            if(host_Gran_neighbour_TVD(1, i) > host_N_cell) then
                print*, "ERROR host_Gran_neighbour_TVD"
                STOP
            end if

            ! if(host_Gran_neighbour_TVD(1, i) < 0) then
            !     print*, "ERROR host_Gran_neighbour_TVD", host_Gran_neighbour_TVD(1, i)
            !     STOP
            ! end if

            ! if(host_Gran_neighbour_TVD(2, i) < 0) then
            !     print*, "ERROR host_Gran_neighbour_TVD", host_Gran_neighbour_TVD(2, i)
            !     STOP
            ! end if
        end do

        call flush(6)

    end subroutine Set_Storage


    subroutine Save_Storage()
        integer :: unit, ierr
        real(8) :: cf

        !3.0  0 градусов
        !3.1  10 градусов
        !3.2  30 градусов


        ! Открываем файл для записи в бинарном формате
        open(newunit=unit, file="FCMHD_1.4_out.bin", form='unformatted', access='stream', &
            action='write', status='replace', iostat=ierr)
        
        if (ierr /= 0) then
            print *, "Error opening file for writing: ", "FCMHD_3.0_out.bin"
            return
        endif
        
        ! Записываем данные в ТОМ ЖЕ порядке, что и при чтении
        write(unit) host_time_all
        write(unit) host_Cell_par 

        cf = 321.0_8
        write(unit) cf

        close(unit)
    end subroutine Save_Storage

    subroutine Fill_data()
        integer :: unit, ierr
        real(8) :: cf

        ! Открываем файл для записи в бинарном формате
        open(3, file = "FCMHD_1.3_out.bin", FORM = 'BINARY', ACTION = "READ")
        
        ! Записываем данные в ТОМ ЖЕ порядке, что и при чтении
        read(3) host_time_all
        read(3) host_Cell_par 

        read(3) cf
        print*, "Proverka (321)  ", cf

        close(3)
        call flush(6)


    end subroutine Fill_data


end module STORAGE