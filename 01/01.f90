! ======================================================================
!     数値流体ゼミ 課題 1
!     波動方程式を数値的に解く.
!
!     ∂u/∂t + c*∂u/∂x = 0 eq:2.6
!     c    = 1
!     CFL  = 0.5
!     Grid = 100 mesh
!     I.C.
!     u(0,x) = 1 (x <= 0.5)
!     u(0,x) = 0 (x >  0.5)
!     B.C.
!     u[0]      = u[1]    : left boundary
!     u[imax+1] = u[imax] : right boundary
!
!     numerical flux 数値流速 : 1 次精度風上差分
!     ~f_i+1/2 = 1/2{(f_i+1 + f_i) - |c|(u_i+1 - u_i)} eq:2.51
!
!     integration : Euler 陽解法
!     u^n+1_i = u^n_i - (∆t/∆x)(~f_i+1/2 - ~f_i-1/2)  eq:2.45
!
!     Code Flow :
!     1. 格子生成
!     2. 初期条件を設定
!     3. 数値流速を計算
!     4. 時間積分
!     5. 境界処理
!     6. 3→6 をループ
!     7. 終了
! ======================================================================

! =====================
! module
! =====================
module myMod
    implicit none
    ! Note: change 13 -> 5 to make everything single precision
    integer , parameter :: p2 = selected_real_kind(13) ! Double Precision
    real(p2), parameter :: zero = 0.0_p2
    real(p2), parameter ::  one = 1.0_p2
    real(p2), parameter :: half = 0.5_p2

    ! セルデータの構造体
    type cell_data
        real(p2) :: xg   ! Cell-center coordinate
        real(p2) :: u    ! Conservative variable
        real(p2) :: u0   ! previous time step value
        real(p2) :: res  ! Residual = f_{j+1/2) - f_{j-1/2)
    end type cell_data

    ! Local variables
    type(cell_data), allocatable :: cell(:)    ! Array of cell-data
    real(p2), allocatable        :: flux(:)    ! Numerical Flux
    real(p2)                     :: xmin, xmax ! 領域の端
    real(p2)                     :: dx, dt     ! 空間・時間刻み幅
    real(p2)                     :: t, tf      ! 現在と最後の時間
    real(p2)                     :: cfl        ! CFL 数と時間刻み幅
    integer                      :: ncells     ! セル数
    integer                      :: nsteps     ! 時間ステップ数
    integer                      :: itime      ! 時間ステップのインデックス
    integer                      :: i, j       ! ループで使う
    integer                      :: step,stepf ! ステップ数と終了ステップ数
    integer                      :: wrtStep    ! 書き出すステップ
end module myMod

! ======================================================================
! main
! ======================================================================
program main
    use myMod
    implicit none

    ! Input params
    call param
    ! grid
    call grid
    ! I.C.
    call init
    ! write Initial profile
    call wrtData

    t = 0.0_p2
    step = 0
    ! loop 時間かステップが終了するまで
    do while ((t < tf).and.(step < stepf))
        step = step + 1
        t    = t + dt
        ! numerical flux
        call numFlux
        ! time integration
        call intg
        ! Boundary Condition
        call bdry
        ! write data
        call wrtData
    end do

    write(6,*)'end! t=',t,'step=',step

end program main
! ======================================================================
! main end
! ======================================================================

! ===============================================
! Input parameters and I.C.
! ===============================================
subroutine param
    use myMod
    ! parameters
    ncells  = 100    ! numbers of cells
    tf      = 1.0_p2 ! Final time
    stepf   = 100    ! Final step
    wrtStep = 5      ! write step
    cfl     = 0.5    ! CFL number
    xmin    = 0.0_p2 ! Left boundary coordinate
    xmax    = 1.0_p2 ! Right boundary coordinate

    ! ゴーストセルも含めて allocate する
    allocate(cell(0:ncells+1)) ! セルデータ
    allocate(flux(1:ncells+1)) ! Numerical Flux

    ! セル幅
    dx = (xmax-xmin)/real(ncells)
    ! 時間刻み幅
    dt = cfl*dx
end subroutine param

! ===============================================
! generate grid
! ===============================================
subroutine grid
    use myMod
    do i = 0, ncells+1
        cell(i)%xg = xmin+real(i-1)*dx
    end do
end subroutine grid


! ===============================================
! Initial Condition
! ===============================================
subroutine init
    use myMod
    do i = 0, ncells+1
        if (cell(i)%xg <= 0.5_p2) then
            cell(i)%u = 1.0_p2
        else
            cell(i)%u = 0.0_p2
        endif
    end do
end subroutine init

! ===============================================
! Numerical Flux
! 1 次精度風上差分
! ~f_i+1/2 = 1/2{(f_i+1 + f_i) - |c|(u_i+1 - u_i)}
! ===============================================
subroutine numFlux
    use myMod
    do i = 1, ncells+1
        flux(i) = cell(i)%u
    end do
end subroutine numFlux

! ===============================================
! Time Integration
! ===============================================
subroutine intg
    use myMod
    do i = 1, ncells
        cell(i)%u = cell(i)%u - (dt/dx)*(flux(i+1)-flux(i))
    end do
end subroutine intg

! ===============================================
! Boundary Condition
! ===============================================
subroutine bdry
    use myMod
    cell(0)%u        = cell(1)%u
    cell(ncells+1)%u = cell(ncells)%u
end subroutine bdry

! ===============================================
! Output data
! ===============================================
subroutine wrtData
    use myMod
    character field*50
    if(mod(step,wrtStep).eq.0) then
        write(field,'(a,i6.6,a)') "./output/field",step,".dat"
        open(10,file=field,form='formatted',status='replace')
        do i = 1,ncells
            write(10,*) (cell(i)%xg,cell(i)%u)
        end do
        close(10)
    end if
end subroutine wrtData
