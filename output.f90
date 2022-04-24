!------------------------------------------------------------------------------!
! G-CTMQC code version 1.0                                                     !
!------------------------------------------------------------------------------!
!                                                                              !
! AUTHOR: Federica Agostini, Institut de Chimie Physique, Paris-Saclay.        !
!                                                                              !
! MODULE: output                                                               !
! Free software: you can redistribute it and/or modify it under the terms of   !
!                the GNU Lesser General Public License as published by the     !
!                Free Software Foundation, either version 3 of the License, or !
!                (at your choice) any later version.                           !
!                                                                              !
! DESCRIPTION: Calls to various subroutines that output results along the      !
! the dynamics.                                                                !
!                                                                              !
! REVISION HISTORY: 11/03/2021 - Initial Version                               !
!                                                                              !
!------------------------------------------------------------------------------!
!> @author
!> Federica Agostini, Institut de Chimie Physique, University Paris-Saclay.
!> @brief Output subroutines that print: electronic populations and coherences
!! as functions of time; electronic coefficients as functions of positions
!! at different time steps; positions, momenta and energies at different time steps.
MODULE output

  USE variables
  USE kinds
  USE analytical_potentials

  IMPLICIT NONE

  CONTAINS

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: General output subroutine                                    !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic populations and coherences in output
  !! and calls additional output subroutines.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] Vcl velocities of the trajectories
  !> @param[in] BOsigma electronic density matrix
  !> @param i,j integer indices
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param index_ij integer index running over the pairs of electronic states
  !> @return Energies, forces and non-adiabatic couplings are stored in the
  !! arrays BOenergy, BOforce, coup.
  SUBROUTINE plot(BOsigma,Rcl,Vcl,time)

    INTEGER         ,INTENT(IN)    :: time
    REAL(KIND=DP)   ,INTENT(IN)    :: Rcl(ntraj,n_dof),Vcl(ntraj,n_dof)
    COMPLEX(KIND=QP),INTENT(IN)    :: BOsigma(ntraj,nstates,nstates)
    INTEGEr                        :: i,j,itraj,index_ij

    IF(time==0 .AND. new_potential) CALL plot_potential

    CALL plot_coefficients(BOsigma,Rcl,time)

    DO itraj=1,ntraj
      CALL compute_energy(BOsigma(itraj,:,:),BOenergy(itraj,:),itraj)
    END DO

    CALL plot_R_P_E(Rcl,Vcl,time)

    IF(time==0) CALL initialize_output

    BO_pop = 0.0_dp
    BO_pop_SH = 0.0_dp
    IF (typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") THEN
       DO i=1,nstates
          DO itraj=1,ntraj
             IF (occ_state(itraj)==i) BO_pop_SH(i) = BO_pop_SH(i) + 1_dp 
          END DO
       END DO
       BO_pop_SH = BO_pop_SH/dble(ntraj)
    END IF
    DO itraj=1,ntraj
       DO i=1,nstates
          BO_pop(i) = BO_pop(i) + real(BOsigma(itraj,i,i),kind=dp)
       END DO
    END DO
  
    BO_pop = BO_pop/dble(ntraj)

    BO_coh   = 0.0_dp
    index_ij = 0

    DO i=1,nstates
       DO j=i+1,nstates
          index_ij=index_ij+1
          DO itraj=1,ntraj
             BO_coh(index_ij) = BO_coh(index_ij) + &
               (real(BOsigma(itraj,i,i),KIND=dp))* &
               (real(BOsigma(itraj,j,j),KIND=dp))
          END DO
       END DO
    END DO
    BO_coh = BO_coh/dble(ntraj)

    WRITE(88,'(f14.4,100f14.8)') DBLE(time)*dt,BO_coh
    WRITE(89,'(f14.4,100f14.8)') DBLE(time)*dt,BO_pop,BO_pop_SH
    
    IF(time==nsteps) CALL finalize_output

  END SUBROUTINE plot

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Output subroutine for electronic coefficients                !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic coeffecients as functions of the
  !! trajectory positions at some selected time steps along the dynamics.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] BOsigma electronic density matrix
  !> @param idx index labelling the output files from 000 to 999
  !> @param filename name of the output file
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param ios control variable for output errors
  !> @return In the directory coeff the files coeff.XXX.dat are created, labelled
  !! from 000 to 999 (those indices label the time steps).
  SUBROUTINE plot_coefficients(BOsigma,Rcl,time)

    INTEGER,INTENT(IN)          :: time
    REAL(KIND=DP),INTENT(IN)    :: Rcl(ntraj,n_dof)
    COMPLEX(KIND=QP),INTENT(IN) :: BOsigma(ntraj,nstates,nstates)
    CHARACTER(LEN=3)            :: idx
    CHARACTER(LEN=400)          :: filename
    INTEGER                     :: itraj,ios

    WRITE(idx,'(i3.3)') time/dump
    filename=TRIM(output_folder)//"/coeff/coeff."//TRIM(idx)//".dat"

    OPEN(128,FILE=TRIM(filename),STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening coefficients file'
    WRITE(128,*) "#Postion, Coefficients: Real part and Imaginary part"

    DO itraj=1,ntraj
      WRITE(128,'(100f14.8)') Rcl(itraj,:), &
        (REAL(BOsigma(itraj,:,:))),(AIMAG(BOsigma(itraj,:,:)))
    END DO    

    CLOSE(128)

  END SUBROUTINE plot_coefficients

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Output subroutine for nuclear properties                     !
  !---------------------------------------------------------------------------!
  !> Subroutine which writes electronic coeffecients as functions of the
  !! trajectory positions at some selected time steps along the dynamics.
  !> @param[in] time time step
  !> @param[in] Rcl positions of the trajectories
  !> @param[in] Vcl velocities of the trajectories
  !> @param idx index labelling the output files from 000 to 999
  !> @param filename name of the output file
  !> @param itraj integer index running over the Ntraj trajectories
  !> @param ios control variable for output errors
  !> @return In the directory trajectories the files RPE.XXX.dat are created,
  !! labelled from 000 to 999 (those indices label the time steps).
  SUBROUTINE plot_R_P_E(Rcl,Vcl,time)

    INTEGER      ,INTENT(IN) :: time
    REAL(KIND=DP),INTENT(IN) :: Rcl(ntraj,n_dof), &
                                Vcl(ntraj,n_dof)
    CHARACTER(LEN=3)         :: idx
    CHARACTER(LEN=400)       :: filename
    INTEGER                  :: itraj,ios

    WRITE(idx,'(i3.3)') time/dump
    FILENAME=TRIM(output_folder)//"/trajectories/RPE."//TRIM(idx)//".dat"
    OPEN(128,file=TRIM(filename),STATUS="replace", &
      FORM="formatted",action="write",iostat=ios)
    IF(ios/=0) PRINT*,'error opening RPE file'
    IF (typ_cal=="TSHFS" .OR. typ_cal=="TSHLZ") THEN
       WRITE(128,*) "#Postions, Momenta, BOsurfaces, Running surface"
    ELSE
      WRITE(128,*) "#Postions, Momenta, TDPES (GI part)"
    END IF

    filename=TRIM(output_folder)//"/trajectories/occ_state."//TRIM(idx)//".dat"
    OPEN(138,FILE=TRIM(filename),STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening occ_state file'

    DO itraj=1,ntraj
       IF (typ_cal=="TSHLZ" .OR. typ_cal=="TSHFS") THEN
          WRITE(128,'(200f14.6)') Rcl(itraj,:), mass(:) * Vcl(itraj,:), &
              BOenergy(itraj,occ_state(itraj)),BOenergy(itraj,:),tdpes(itraj)
          WRITE(138,'(i5)') occ_state(itraj),LZ_hop(itraj)
       ELSE
          WRITE(128,'(200f14.6)') Rcl(itraj,:),mass(:) * Vcl(itraj,:),  &
                 tdpes(itraj),BOenergy(itraj,:),coup(itraj,1,2,:)
       END IF
    END DO

    CLOSE(128)
    CLOSE(138)

  END SUBROUTINE plot_R_P_E

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Calculation of the potential energy in Ehrenfest and CT-MQC  !
  !---------------------------------------------------------------------------!
  !> The subroutine computes the expectation value of the electronic Hamiltonian on
  !! the time-dependent electronic wavefunction, yielding the gauge-invariant
  !! part of the TDPES in CT-MQC or the mean Ehrenfest potential.
  !> @param[in] trajlabel label of the trajectory
  !> @param[in] my_rho electronic density matrix
  !> @param[in] e_BO adiabatic or spin-(a)diabatic energy
  !> @param i integer index
  !> @return The value of the TDPES is returned, where "TDPES" means either
  !! the gauge-invariant part of the TDPES in CT-MQC or the mean Ehrenfest
  !! potential.
  SUBROUTINE compute_energy(my_rho,e_BO,trajlabel)

    INTEGER         ,INTENT(IN) :: trajlabel
    COMPLEX(KIND=QP),INTENT(IN) :: my_rho(nstates,nstates)
    REAL(KIND=DP)   ,INTENT(IN) :: e_BO(nstates)
    INTEGER                     :: i

    tdpes(trajlabel)=0.0_dp

    DO i=1,nstates
      tdpes(trajlabel)=tdpes(trajlabel)+ &
        REAL(my_rho(i,i),KIND=dp)*e_BO(i)
    END DO

  END SUBROUTINE compute_energy

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Initialization of the output                                 !
  !---------------------------------------------------------------------------!
  !> The files where electronic populations and coherences are written are
  !! opened in this subroutine.
  !> @param ios control variable for output errors
  SUBROUTINE initialize_output

    INTEGER :: ios

    OPEN(88,FILE=TRIM(output_folder)//"/BO_coherences.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_coherences.dat'
    WRITE(88,*) "#Time, Coherences"

    open(89,FILE=TRIM(output_folder)//"/BO_population.dat",STATUS="replace", &
      FORM="formatted",ACTION="write",IOSTAT=ios)
    IF(ios/=0) PRINT*,'error opening BO_population.dat'
    WRITE(89,*) "#Time, Populations"

  END SUBROUTINE initialize_output

  !---------------------------------------------------------------------------!
  ! DESCRIPTION: Finalization of the output                                   !
  !---------------------------------------------------------------------------!
  !> The files where electronic populations and coherences are written are
  !! closed in this subroutine.
  SUBROUTINE finalize_output

    close(89)
    close(88)

  END SUBROUTINE finalize_output

END module output
