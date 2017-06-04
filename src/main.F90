program main

  use constants
  use finalize,          only: openmc_finalize
  use global
  use initialize,        only: openmc_init
  use message_passing
  use particle_restart,  only: run_particle_restart
  use plot,              only: run_plot
  use simulation,        only: run_simulation
  use volume_calc,       only: run_volume_calculations

  use physics,           only: sample_cxs_target_velocity

  implicit none

  integer :: i
  real(8) :: v_target(3)

  ! Initialize run -- when run with MPI, pass communicator
#ifdef MPI
  call openmc_init(MPI_COMM_WORLD)
#else
  call openmc_init()
#endif

  !! start problem based on mode
  !select case (run_mode)
  !case (MODE_FIXEDSOURCE, MODE_EIGENVALUE)
  !  call run_simulation()
  !case (MODE_PLOTTING)
  !  call run_plot()
  !case (MODE_PARTICLE)
  !  if (master) call run_particle_restart()
  !case (MODE_VOLUME)
  !  call run_volume_calculations()
  !end select

  open(unit=2, file="vel.txt")

  do i=1, 100000
    call sample_cxs_target_velocity(nuclides(1), v_target, ZERO, &
                                    [ONE, ZERO, ZERO], K_BOLTZMANN * 300.0_8)
    write(2, "(e10.3)") &
        sqrt(TWO * dot_product(v_target, v_target) / MASS_NEUTRON_EV) * C_LIGHT
  end do

  close(2)

  !! finalize run
  !call openmc_finalize()

end program main
