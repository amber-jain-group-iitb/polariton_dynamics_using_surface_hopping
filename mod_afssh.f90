Module mod_afssh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28,ev_2_wave=8065.73
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J,omg_max,w_c
integer,parameter :: g=500

!! Potential
integer nquant
real*8 g_coup,epsilon
real*8 V_exothermicity,Vc,omg_B,gamma_B,temperature
real*8 beta,gamma_D,lambda_B,V_reorg,V_barrier,V_coup
real*8 omg_c,omg_scaled
real*8 s01,s02,x_cr
real*8,allocatable :: mass(:),omg(:),ck(:)

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)
complex*16,allocatable :: rho(:,:,:),rho1(:),rho2(:),rho3(:),rho4(:),rho5(:),rho6(:),rho7(:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:,:),v(:,:),acc(:,:)
real*8,allocatable :: x_old(:,:),v_old(:,:),acc_old(:,:),x_hop(:,:)
real*8 tim_hop
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
complex*16,allocatable :: delr(:,:,:,:),delp(:,:,:,:),delacc(:,:,:,:)
complex*16,allocatable :: delr_old(:,:,:,:),delp_old(:,:,:,:)

!! Quantized vibration
integer ncl_site
integer nb_vib,n_dvr
real*8,allocatable ::si_sho(:,:,:),sho_overlap(:,:,:,:),q_exp(:,:,:)
real*8,allocatable ::qsq_exp(:,:,:)
real*8,allocatable ::en_sho(:,:),fc_init(:)
real*8,allocatable ::Hamil_diab_0(:,:)

!! Quantum
integer state,nbasis,state_tentative
integer state_old
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_site(:,:),Hamil_diab(:,:),delH_dels(:,:,:,:),delH_dels_ad(:,:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:,:),force_old(:,:,:,:),delF(:,:,:,:)
complex*16,allocatable :: ci(:),ci_old(:),sigma(:,:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:),W_overlap(:,:),hop_prob_net(:)
real*8 x_grid(g),sho_hamil_1(g,g),sho_hamil_2(g,g),sho_hamil_0(g,g)
real*8 expec_q11(g,g),expec_q22(g,g),expec_q12(g,g),expec_q21(g,g)
real*8 x_grid_max,x_grid_min,dx,ens1(g),vect1(g,g),ens2(g),vect2(g,g),ens0(g),vect0(g,g)
real*8 expec_sqr_q11(g,g),expec_sqr_q22(g,g),expec_sqr_q12(g,g),expec_sqr_q21(g,g)
real*8 expec11(g,g),expec22(g,g),expec12(g,g),expec21(g,g)

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

!!variable for photons
integer::n_photon
real*8::omg_photon,coup_intra,coup_photon,coup_el
real*8,allocatable::site_en(:)

!!explicit bath variables
integer::nbath,nsite,DSE
real*8,allocatable::omg_bath(:,:),ci_bath(:,:)
real*8::bath_freq_min,bath_freq_max,del_bath_freq,mass_b

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2),seed_out
  real*8 rnd,c_0,c_e,kt

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) flag_ortho
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) n_dvr
  read(10,*) nb_vib
  read(10,*) V_exothermicity
  read(10,*) Vc
  read(10,*) omg_B
  read(10,*) lambda_B
  read(10,*) gamma_B
  read(10,*) temperature
  read(10,*) nbath
  read(10,*) DSE
  read(10,*) coup_photon
  read(10,*) bath_freq_min
  read(10,*) bath_freq_max
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  V_coup=-87.7
  n_photon=1
  nsite=nquant
  nquant=nquant+n_photon
  ncl_site=nclass
  nbasis=nquant
  bath_freq_min=bath_freq_min*2*pi*clight
  bath_freq_max=bath_freq_max*2*pi*clight

  coup_photon=coup_photon*wave_to_J
  energy_cutoff=energy_cutoff*wave_to_J
!  temperature=temperature*wave_to_J/kb
  kt=kb*temperature
  Vc=Vc*wave_to_J
  V_exothermicity=V_exothermicity*wave_to_J
  omg_B=omg_B*2*pi*clight
  gamma_B=gamma_B*2*pi*clight
  lambda_B=lambda_B*wave_to_J

  nsteps=nint(total_time/dtc)+1
  beta=1.d0/(kb*temperature)
  !lambda_D=V_reorg/4.d0

  !-----------------------------------------------------------------  
  i=(nsteps)/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(rho(nquant,nquant,i),rho1(i),rho2(i),rho3(i),rho4(i),rho5(i),rho6(i),rho7(i))
  allocate(x(nbath,nsite),v(nbath,nsite),acc(nbath,nsite))
  allocate(x_old(nbath,nsite),v_old(nbath,nsite),acc_old(nbath,nsite),x_hop(nbath,nsite))
  allocate(mass(nsite),omg(nsite),ck(nsite))
  allocate(delr(nquant,nquant,nbath,nsite),delp(nquant,nquant,nbath,nsite),delacc(nquant,nquant,nbath,nsite))
  allocate(delr_old(nquant,nquant,nbath,nsite),delp_old(nquant,nquant,nbath,nsite))
  allocate(si_adiab(nbasis,nquant),ci(nquant),V_k(nquant),V_k_old(nquant),sigma(nquant,nquant))
  allocate(Hamil_site(nbasis,nbasis),Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nbath,nsite),delH_dels_ad(nquant,nquant,nbath,nsite))
  allocate(pot(nquant,nquant),force(nquant,nquant,nbath,nsite),force_old(nquant,nquant,nbath,nsite),delf(nquant,nquant,nbath,nsite))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nbath,nsite),vdotd(nquant,nquant),hop_prob(nquant),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant),si_adiab_prev(nbasis,nquant))
  allocate(si_sho(n_dvr,nb_vib,2))
  allocate(sho_overlap(nb_vib,nb_vib,2,2))
  allocate(q_exp(nb_vib,nb_vib,2),qsq_exp(nb_vib,nb_vib,2))
  allocate(en_sho(nb_vib,2),fc_init(nb_vib))
  allocate(Hamil_diab_0(nbasis,nbasis))
  allocate(site_en(nsite),omg_bath(nbath,nsite),ci_bath(nbath,nsite))

open(20,file='ifolder.inp')
read(20,*)seed_out
close(20)

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)*seed_out
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)*seed_out
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        !call init_cond
        call random_number(rnd)
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0

end subroutine setup
!---------------------------------------------------------- 
subroutine main
  implicit none
  integer i,j,k,n
  real*8 Vc_save
integer::t1,t2
INTEGER :: count_rate

  call files(0)

!  call cpu_time(t1)
call system_clock(t1,count_rate)

  call setup_parameters
  call initialize_averages
  call bath_coupling

 Vc_save=V_coup

!call check_acceleration
do i=1,N_traj
    traj_num=i
    !V_coup=0.001d0
    !call setup_parameters
    call init_cond
    !call evolve(1,5000)
    !V_coup=Vc_save;call setup_parameters
    call evolve(1,nsteps)
    call average_end
  enddo
call write_average

 call system_clock(t2,count_rate);tim_tot=(tim_tot+t2-t1)/count_rate
 call files(1)

end subroutine main

!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_qm")
    open(13,file="output_hop")
    open(14,file="output_overlap")
    open(15,file="output_dec")

    open(100,file="pop1.out")
    open(101,file="cnts.out")
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    close(10);close(11);close(12);close(13);close(14);close(15)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  cnt_frust=0.d0
  cnt_collapse=0.d0
  cnt_init=0.d0
  cnt_term=0.d0
  pop=0.d0
  rho=0.d0
  rho1=0.d0
  rho2=0.d0
  rho3=0.d0
  rho4=0.d0
  rho5=0.d0
  rho6=0.d0
  pop_surf=0.d0
  pop_amp=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  
subroutine init_cond
  implicit none
  integer i,j
  real*8 sig_x,sig_p,rnd,ak,su
  real*8 energy_0,part_fnc
complex*16::ci_adiab(nquant,1),rho_tmp(nquant,nquant),ci_diab(nquant,nquant)

!  do i=1,nclass
    !ak=2/(hbar*omg(i))*dtanh(beta*hbar*omg(i)/2.d0) !! Wigner
!    ak=beta    !! Classical
!    sig_x=1.d0/dsqrt(ak*mass(i)*omg(i)**2)
!    sig_p=dsqrt(mass(i)/ak)
!    call gaussian_random_number(rnd)
!    x(i)=rnd*sig_x
!    call gaussian_random_number(rnd)
!    v(i)=(1.d0/mass(i)*(rnd*sig_p))
!  enddo
  state=1

do i=1,nsite
   do j=1,nbath
        ak=beta
        sig_x=1.d0/dsqrt(ak*mass_b*omg_bath(j,i)**2)
        sig_p=dsqrt(mass_b/ak)
        call gaussian_random_number(rnd)
        x(j,i)=rnd*sig_x
!        x(j,i)=x(j,i)+ci_bath(j,i)/(mass_b*omg_bath(j,i)**2)
        call gaussian_random_number(rnd)
        v(j,i)=1.d0/mass_b*rnd*sig_p
   enddo
enddo

!call equilibrate

  call evaluate_variables(0)
  call evaluate_variables(1)

ci(:)=si_adiab(state,:)

 call random_number(rnd)
  su=0.d0
  do i=1,nquant
    su=su+cdabs(ci(i))**2
    if(rnd<su) then
      state=i
      exit
    endif
  enddo

!write(*,*) state
  delr=0.d0
  delp=0.d0

  ihop=1
  iaverage=1
  iterminate=0
  flag_terminate=0

  curr_time=0.d0
  call evaluate_variables(0)
  call evaluate_variables(1)
  call compute_mat_diab

  !! to compute the standard deviation of the energy of the trajectory
  en_avg=0.d0;ensq_avg=0.d0
  nst_av=0

end subroutine init_cond

!-----------------------------------------------------------------  
subroutine equilibrate
implicit none
integer::i,j
real*8::dte,energy1

dte=1.d-15
do i=1,500000
x=x+v*dte+0.5*acc*dte*dte
v=v+acc*dte*0.5
call acc_eq
v=v+acc*dte*0.5

!energy1=0.d0
!do j=1,nbath
!energy1=energy1+0.5d0*mass_b*(omg_bath(j,1)**2*(x(j,1)-ci_bath(j,1)/(mass_b*omg_bath(j,1)**2))**2+omg_bath(j,2)**2*(x(j,2)-ci_bath(j,1)/(mass_b*omg_bath(j,2)**2))**2)
!energy1=energy1+0.5d0*mass_b*(v(j,1)**2+v(j,2)**2)
!enddo
!write(50,*)i,energy1/wave_to_J
enddo

!stop

end subroutine
!-----------------------------------------------------------------  
subroutine acc_eq
implicit none
integer::i,j

do j=1,nbath
do i=1,nsite
acc(j,i)=-omg_bath(j,i)**2*(x(j,i)-ci_bath(j,i)/(mass_b*omg_bath(j,i)**2))
enddo
enddo

end subroutine
!-----------------------------------------------------------------  
 subroutine evolve(start,nsteps)
  implicit none
  integer,intent(in) :: start,nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  call write_output(1,1)
  iterm=0
  do i=start,nsteps
    call write_output(i,0)
    call average(i)
    call save_old_state
    call evolve_classical(dtc)
    i_do_something=0


    if(ifriction==0) then
      !! Do till energy is conserved
      !! ifriction==1 --> Langevin equation (non-energy conserving calculations)
      do while(dabs(energy-energy_old)>energy_cutoff.and.i>1) 
        i_do_something=i_do_something+1
        call do_something(i_do_something)
      enddo
    endif
    if(i_do_something==1) then
      !! if_do_something==1 refers to energy conservation by changing adiabat.
      !! This is treated as if a hop occured - hence reset moments
      delr=0.d0
      delp=0.d0
    endif

!write(200,*)curr_time*1.d15,i_do_something
    if(i_do_something.ne.1)call evolve_quantum_small_dtq
    !call evolve_quantum_small_dtq
    if(i_do_something.ne.1.and.ihop==1)call hop
    if(i_do_something.ne.1.and.icollapse==1)call collapse(dtc,iflag_coll)
    if(flag_terminate==1) call traj_terminate(iterm)
      if(iterm==1)exit

    curr_time=curr_time+dtc
  enddo
  call write_output(1,1)

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine do_something(i_do_something)
  !! subroutine to maintain energy conservation
  implicit none
  integer,intent(in)::i_do_something
  real*8 acc_tent(nclass),dt,dtm(3)
  integer i,nstep

  if(i_do_something==1) then
    !! On first pass, check if hop happens; if yes, check if evolution using the acceleration of hopped surface conservses energy
    !! Useful for very sharp crossings
    call evolve_quantum_small_dtq
    if(flag_hop==1) then
      state=state_tentative
      call evaluate_variables(0)
      v=v_old+0.5*(acc_old+acc)*dtc
      call evaluate_variables(1)
    endif
  else
    !! If the first pass did not work, reduce time-steps untill energy is conserved
    dtm=1.d0 !! some large value
    dtm(1)=0.1d0/maxval(vdotd)
    dtm(2)=0.5*dtc*dsqrt(energy_cutoff/dabs(energy-energy_old))
    dtm(3)=dtc

    dt=minval(dtm)

    dt=dt/real(i_do_something)
    nstep=nint(dtc/dt)
    dt=dtc/real(nstep)
    call revert_state
    do i=1,nstep
      call evolve_classical(dt)
    enddo
  endif


end subroutine do_something
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer j,i1,j1,k,kp
  complex*16 ci_diab(nquant,nquant),rho_ad(nquant,nquant),ci_adiab(nquant,1)
  real*8 r_avg,U(nquant,nquant),U_exc(nquant,nquant)
  integer if_reactant
  real*8 t1,t2

  !call cpu_time(t1)

  if(iwrite==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
    nst_av=nst_av+1
  endif

  if(iaverage==1.and.(mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif

    U=si_adiab
    rho_ad=0.d0
    rho_ad(state,state)=1.d0
    do i1=1,nquant
      do j1=1,nquant
        if(i1.ne.j1) rho_ad(i1,j1)=ci(i1)*dconjg(ci(j1))
      enddo
    enddo
    rho(:,:,j)=rho(:,:,j)+matmul(U,matmul(rho_ad,transpose(U)))

!ci_adiab(:,1)=ci(:)
!rho_ad=matmul(conjg(ci_adiab),transpose(ci_adiab))
!ci_diab=matmul(si_adiab,matmul(rho_ad,transpose(si_adiab)))

!do i1=1,nbasis
!        pop(i1,j)=pop(i1,j)+cdabs(ci_diab(i1,i1))
!enddo
    
    !pop(:,j)=pop(:,j)+si_adiab(:,state)**2
    !pop_surf(:,j)=pop_surf(:,j)+si_adiab(:,state)**2
    !ci_diab=matmul(si_adiab,ci)
    !pop_amp(:,j)=pop_amp(:,j)+cdabs(ci_diab)**2
    !do j1=2,nquant
    !  do i1=1,j1-1
    !    pop(:,j)=pop(:,j)+2*real(ci(i1)*dconjg(ci(j1)))*si_adiab(:,i1)*si_adiab(:,j1)
    !  enddo
    !enddo
  endif

  !call cpu_time(t2)
  !tim_coll=tim_coll+t2-t1

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

end subroutine average_end
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  acc_old=acc
  ci_old=ci
  state_old=state
  !ci2_old=ci2
  si_adiab_prev=si_adiab
  V_k_old=V_k
  force_old=force
  energy_old=energy
  delr_old=delr
  delp_old=delp

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  state=state_old
  ci=ci_old
  delr=delr_old
  delp=delp_old
  force=force_old
  !ci2=ci2_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum_small_dtq
  implicit none
  integer i,nstep_qm
  real*8 dtq,dtq1,dtq2
  real*8 V_k_hold(nquant),dVk_dt(nquant)
  complex*16 ci_prev(nquant),dci_dt(nquant)

  call compute_vdotd
  dVk_dt=(V_k-V_k_old)/dtc
  if(icollapse==1) then
    call compute_delH_dels_ad
  endif

  dtq1=0.02/maxval(vdotd)
  dtq2=0.02*hbar/maxval(V_k-sum(V_k)/real(nquant))
  dtq=dtq1
  if(dtq>dtq2)dtq=dtq2

  if(dtq>dtc)dtq=dtc
  nstep_qm=nint(dtc/dtq)
  dtq=dtc/real(nstep_qm)
  hop_prob=0.d0
  hop_prob_net=0.d0
  V_k_hold=V_k
  V_k=V_k_old
  call compute_mat_adiab


  flag_hop=0
  do i=1,nstep_qm
    call compute_hop_prob(dtq)
    if(flag_hop==0)call check_hop(i*dtq)
    call rk4(ci,dtq,dVk_dt)
  !  if(icollapse==1)call rk4_decoherence(dtq)
  enddo
  if(icollapse==1)call vv_decoherence(dtc)

  !if(icollapse==1) then
  !  call verlet_decoherence(dtc,W_overlap,V_k_old,dvk_dt)
  !endif

  do i=1,nquant
    if(hop_prob_net(i)<0.d0)hop_prob_net=0.d0
    hop_prob_net(i)=1.d0-dexp(-hop_prob_net(i))
  enddo

end subroutine evolve_quantum_small_dtq
!-----------------------------------------------------------------  

subroutine compute_hop_prob(dtq)
  implicit none
  real*8,intent(in)::dtq
  integer i
  real*8 pr

  do i=1,nquant
    if(i.ne.state) then
      pr=-2*real(ci(i)*dconjg(ci(state)))*vdotd(i,state)
      pr=pr*dtq/cdabs(ci(state))**2
      if(pr<0.d0)pr=0.d0     !!!! CAUTION AMBER CHECK !!!!
      hop_prob(i)=pr
      hop_prob_net(i)=hop_prob_net(i)+pr
    endif
  enddo

end subroutine compute_hop_prob
!-----------------------------------------------------------------  

subroutine check_hop(tim)
  implicit none
  real*8,intent(in)::tim
  integer i
  real*8 rnd,pr

  call random_number(rnd)
  pr=0.d0
  flag_hop=0
  do i=1,nquant
    if(i.ne.state) then
      pr=pr+hop_prob(i)
      if(rnd<pr) then
        state_tentative=i
        flag_hop=1
        exit
      endif
    endif
  enddo

end subroutine check_hop
!-----------------------------------------------------------------  

subroutine rk4(ci,dtq,dVk_dt)
  implicit none
  complex*16,intent(inout)::ci(nquant)
  real*8,intent(in) :: dtq,dVk_dt(nquant)
  complex*16,dimension(1:nquant):: k1,k2,k3,k4

  k1=matmul(mat_adiab,ci)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k2=matmul(mat_adiab,ci+0.5*dtq*k1)
  k3=matmul(mat_adiab,ci+0.5*dtq*k2)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k4=matmul(mat_adiab,ci+dtq*k3)

  ci=ci+dtq/6.d0*(k1+2*k2+2*k3+k4)


!write(60,'(20es)')curr_time*1.d15,cdabs(ci)

end subroutine rk4
!-----------------------------------------------------------------  

subroutine vv_decoherence(dtc)
  implicit none
  integer i,j
  real*8,intent(in)::dtc
  real*8 delacc_old(nquant,nquant,nbath,nsite)

  delr=delr+delp/mass(1)*dtc+0.5*delacc*dtc**2/mass(1)
  delacc_old=delacc
  call compute_delacc
  delp=delp+0.5*(delacc+delacc_old)*dtc

  do i=1,nbath
  do j=1,nsite
    delr(:,:,i,j)=matmul(W_overlap,matmul(delr(:,:,i,j),W_overlap))
    delp(:,:,i,j)=matmul(W_overlap,matmul(delp(:,:,i,j),W_overlap))
  enddo
  enddo

  do i=1,nbath
  do j=1,nsite
    delr(:,:,i,j)=delr(:,:,i,j)-delr(state,state,i,j)
    delp(:,:,i,j)=delp(:,:,i,j)-delp(state,state,i,j)
  enddo
  enddo


end subroutine vv_decoherence
!-----------------------------------------------------------------

subroutine compute_delacc
  implicit none
  integer i,j

  do i=1,nbath
  do j=1,nsite
    delacc(:,:,i,j)=0.5*anti_commute(delF(:,:,i,j),sigma,0)
  enddo
  enddo

end subroutine compute_delacc
!-----------------------------------------------------------------  

subroutine rk4_decoherence(dtq)
  implicit none
  real*8,intent(in)::dtq
  complex*16,dimension(2,nquant,nquant,nbath,nsite):: kd1,kd2,kd3,kd4,vec

  vec(1,:,:,:,:)=delr
  vec(2,:,:,:,:)=delp

  call compute_T_jk(kd1,vec)
  call compute_T_jk(kd2,vec+0.5*dtq*kd1)
  call compute_T_jk(kd3,vec+0.5*dtq*kd2)
  call compute_T_jk(kd4,vec+dtq*kd3)

  vec=vec+dtq/6.d0*(kd1+2*kd2+2*kd3+kd4)
  delr=vec(1,:,:,:,:)
  delp=vec(2,:,:,:,:)

end subroutine rk4_decoherence
!-----------------------------------------------------------------  

subroutine compute_T_jk(T_jk,vec)
  implicit none
  complex*16,intent(in):: vec(2,nquant,nquant,nbath,nsite)
  complex*16,intent(out):: T_jk(2,nquant,nquant,nbath,nsite)
!  complex*16 delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16 Tr(nquant,nquant,nbath,nsite)
  complex*16 tmp1(nbath,nsite),tmp2(nbath,nsite)
  integer i
  real*8 t1,t2,t11,t12,t21,t22

  !call cpu_time(t1)

  delr=vec(1,:,:,:,:)
  delp=vec(2,:,:,:,:)

!  call cpu_time(t11)
  call compute_T_jk_R(Tr,delr,delp)
!  call cpu_time(t21);tim_check=tim_check+(t21-t11)

  T_jk(1,:,:,:,:)=Tr
!  call cpu_time(t12)
  call compute_T_jk_P(Tr,delr,delp)
!  call cpu_time(t22);tim_check2=tim_check2+(t22-t12)

  T_jk(2,:,:,:,:)=Tr

  tmp1=T_jk(1,state,state,:,:)
  tmp2=T_jk(2,state,state,:,:)

  do i=1,nquant
    T_jk(1,i,i,:,:)=T_jk(1,i,i,:,:)-tmp1
    T_jk(2,i,i,:,:)=T_jk(2,i,i,:,:)-tmp2
  enddo

  !call cpu_time(t2)
  !tim_T_jk=tim_T_jk+t2-t1

end subroutine compute_T_jk
!-----------------------------------------------------------------  

subroutine compute_T_jk_R(T_jk,delr,delp)
  !! Eq. 14 of JCP 137, 22A513
  implicit none
  complex*16,intent(in) ::  delr(nquant,nquant,nbath,nsite),delp(nquant,nquant,nbath,nsite)
  complex*16,intent(out) :: T_jk(nquant,nquant,nbath,nsite)
  integer i1,j1

  do i1=1,nbath
  do j1=1,nsite
      T_jk(:,:,i1,j1)=-iota/hbar*commute(pot,delr(:,:,i1,j1),0)+delp(:,:,i1,j1)/mass(i1)
      T_jk(:,:,i1,j1)=T_jk(:,:,i1,j1)-commute(vdotd,delr(:,:,i1,j1),0)
  enddo
  enddo

end subroutine compute_T_jk_R
!-----------------------------------------------------------------  

subroutine compute_T_jk_P(T_jk,delr,delp)
  !! Eq. 16 of JCP 137, 22A513
  implicit none
  complex*16,intent(in) ::  delr(nquant,nquant,nbath,nsite),delp(nquant,nquant,nbath,nsite)
  complex*16,intent(out) :: T_jk(nquant,nquant,nbath,nsite)
  real*8 delF(nquant,nquant,nbath,nsite)
  integer i1,j1

  delF=force
  do i1=1,nquant
    delF(i1,i1,:,:)=delF(i1,i1,:,:)-force(state,state,:,:)
  enddo

  do i1=1,nsite
     do j1=1,nbath
      T_jk(:,:,j1,i1)=-iota/hbar*commute(pot,delp(:,:,j1,i1),0)
      T_jk(:,:,j1,i1)=T_jk(:,:,j1,i1)+0.5*anti_commute(delF(:,:,j1,i1),sigma,0)
      T_jk(:,:,j1,i1)=T_jk(:,:,j1,i1)-commute(vdotd,delp(:,:,j1,i1),0)
  enddo
  enddo

end subroutine compute_T_jk_P
!-----------------------------------------------------------------  

!subroutine verlet_decoherence(dt,W_mat,V_k0,dvk_dt)
! implicit none
! real*8,intent(in):: dt,W_mat(nquant,nquant),V_k0(nquant),dvk_dt(nquant)
! real*8 acc_dec(nquant,nclass),delf(nquant,nclass),temp(nclass)
! complex*16 temp_delr(nquant,nclass),temp_delp(nquant,nclass)
! !complex*16 ci_diab(nquant)
! integer i,j,k
!
! delF=force_old
! temp=delF(state,:)
! do i=1,nquant
!   delF(i,:)=delF(i,:)-temp
!   acc_dec(i,:)=delF(i,:)*cdabs(ci_old(i))**2/mass
! enddo
!
! do i=1,nquant
!   delr(i,:)=delr(i,:)+delp(i,:)/mass*dt+0.5*acc_dec(i,:)*dt**2
!   delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dt
! enddo
!
! !ci_diab=cdexp(iota*V_k0*dt/hbar)*cdexp(0.5*iota*dvk_dt*dt**2/hbar)*ci
! !ci_diab=matmul_lap(W_mat,ci_diab)
! delF=0.d0
! do j=1,nquant
!   do k=1,nquant
!     delF(j,:)=delF(j,:)+dabs(W_mat(j,k)**2)*(force(k,:)-force(state,:))
!   enddo
! enddo
! !temp=delF(state,:)
! do i=1,nquant
! !  delF(i,:)=delF(i,:)-temp
! !  !acc_dec(i,:)=delF(i,:)*cdabs(ci_diab(i))**2/mass
!   acc_dec(i,:)=delF(i,:)*cdabs(ci_old(i))**2/mass
! enddo
!
! do i=1,nquant
!   delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dt
! enddo
!
! temp_delr=0.d0;temp_delp=0.d0
! do j=1,nquant
!   do k=1,nquant
!     temp_delr(j,:)=temp_delr(j,:)+dabs(W_mat(k,j)**2)*delr(k,:)
!     temp_delp(j,:)=temp_delp(j,:)+dabs(W_mat(k,j)**2)*delp(k,:)
!   enddo
! enddo
! delr=temp_delr
! delp=temp_delp
!
! !do i=1,nclass
! !  delr(:,i)=delr(:,i)-delr(state,i)
! !  delp(:,i)=delp(:,i)-delp(state,i)
! !enddo
!
!end subroutine verlet_decoherence
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  integer i
  real*8,intent(in) :: dt
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass),acc_sav(nclass)
  real*8 t1,t2

  !call cpu_time(t1)

  if(ifriction==0) then
    !! Step 1
    x=x+v*dt+0.5*acc*dt*dt
    v=v+0.5*acc*dt
    acc_old=acc
    call evaluate_variables(0)
    v=v+0.5*dt*acc
    call evaluate_variables(1)
  endif

  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
!     call stochastic_force(delta_r,delta_v,dt)
!     x=x+c1*dt*v+c2*dt*dt*acc+delta_r
!     acc_old=acc
!     call evaluate_variables(0)
!     v=c0*v+(c1-c2)*dt*acc_old+c2*dt*acc+delta_v
!     call evaluate_variables(1)
  endif

  !call cpu_time(t2);tim_ev_cl=tim_ev_cl+t2-t1

end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine deriv_xv(vec,acc,kk)
  implicit none
  real*8,intent(in)::vec(2*nclass),acc(nclass)
  real*8,intent(out)::kk(2*nclass)

  kk(1:nclass)=vec(nclass+1:2*nclass)
  kk(nclass+1:2*nclass)=acc

end subroutine deriv_xv
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine compute_mat_diab
  implicit none
  integer i,j
  real*8 t1,t2

  !call cpu_time(t1)

  mat=0.d0
  do i=1,nbasis
    do j=1,nbasis
      mat(i,j)=-iota/hbar*sum(si_adiab(i,:)*si_adiab(j,:)*V_k(1:nquant))
    enddo
  enddo

  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1

end subroutine compute_mat_diab
!-----------------------------------------------------------------  

subroutine compute_mat_adiab
  implicit none
  integer i,j
  real*8 t1,t2
  real*8 V_avg
  
  !call cpu_time(t1)

  mat_adiab=-vdotd
  V_avg=sum(V_k)/real(nquant)
  do i=1,nquant
    mat_adiab(i,i)=mat_adiab(i,i)-iota/hbar*(V_k(i)-V_avg)
  enddo
      
  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1
  
end subroutine compute_mat_adiab
!-----------------------------------------------------------------  

subroutine hop
  implicit none
  integer ifrust

  if(flag_hop==1) then
    call velocity_adjust(state_tentative,ifrust)
  endif

end subroutine hop
!-----------------------------------------------------------------  

subroutine velocity_adjust(state_tentative,ifrust)
  implicit none
  integer,intent(in)::state_tentative
  integer,intent(out)::ifrust
  real*8 gij,gama,aa,bb,cc,discr,dp(nbath,nsite),vd,f1,f2
  integer i,j,k,kp

  k=state;kp=state_tentative
  cc=V_k(state)-V_k(state_tentative)

  call compute_dij_2state(x,k,kp,dp)
  dp=dp/dsqrt(sum(dp*dp))

  aa=0.d0
  bb=0.d0
  do i=1,nsite
  do j=1,nbath

    aa=aa+0.5/mass_b*(dp(j,i)*dp(j,i))
    bb=bb+(v(j,i)*dp(j,i))

  enddo
  enddo


  discr=bb**2+4*aa*cc

  if(discr<0.d0) then
    ifrust=1
    cnt_frust=cnt_frust+1.d0
    if(flag_frust==0)then
      gama=0.d0
      call compute_delH_dels_ad
      f1=sum(force(k,k,:,:)*dp)
      f2=sum(force(kp,kp,:,:)*dp)
      vd=sum(v*dp)
      !! reverse velocity based on Truhlar's ideas
      !if(f1*f2<0.d0.and.vd*f2<0.d0) then
      if(vd*f2<0.d0) then
        gama=bb/aa
      endif
    endif
    if(flag_frust>0)gama=0.d0
  else
    ifrust=0
    if(bb>=0.d0) gama=(bb-dsqrt(discr))/(2*aa)
    if(bb<0.d0)  gama=(bb+dsqrt(discr))/(2*aa)
    state=state_tentative
    delr=0.d0
    delp=0.d0
  endif


  do i=1,nsite
  do j=1,nbath
    v(j,i)=v(j,i)-gama*dp(j,i)/mass_b
  enddo
  enddo



!write(20,*)curr_time*1.d15,dp/dsqrt(sum(dp*dp)),x(1),ifrust
!write(21,*)curr_time*1.d15,k,kp,gama

  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine velocity_adjust
!-----------------------------------------------------------------  

subroutine reverse_velocity
  implicit none
  

end subroutine reverse_velocity
!-----------------------------------------------------------------  

subroutine collapse(dt,iflag_coll)
  implicit none
  real*8,intent(in) :: dt
  integer,intent(out) :: iflag_coll
  real*8 rnd,gama_collapse,gama_reset
  complex*16 su1
  integer n,i,j

  i=state

  if(icollapse==1) then

    iflag_coll=0
    do n=1,nquant
      if(n.ne.state) then
        gama_reset=0.d0
        gama_reset=gama_reset+sum((force(n,n,:,:)-force(i,i,:,:))*dble(delr(n,n,:,:)))/(2*hbar)
        su1=0.d0
        su1=su1+sum(force(i,n,:,:)*delr(n,n,:,:))
        gama_collapse=gama_reset-2/hbar*cdabs(su1)
        gama_collapse=gama_collapse*dt
        gama_reset=-gama_reset*dt
        call random_number(rnd)

        if(rnd<gama_collapse) then
          iflag_coll=1
          cnt_collapse=cnt_collapse+1
          if(icollapse==1) then
            !do j=1,nquant
            !  if(j.ne.n) ci(j)=ci(j)/dsqrt(1-cdabs(ci(n)**2))
            !enddo
            !! Erratum: Landry, Subotnik JCP 137, 229901 (2012)
            ci(i)=ci(i)/cdabs(ci(i))*dsqrt(cdabs(ci(i))**2+cdabs(ci(n))**2)
            ci(n)=0.d0
!write(61,'(19es)')curr_time*1.d15,cdabs(ci)
          endif
        endif
        if(rnd<gama_collapse.or.rnd<gama_reset) then
          if(icollapse==1) then
            do j=1,nquant
              delr(j,n,:,:)=0.d0;delr(n,j,:,:)=0.d0
              delp(j,n,:,:)=0.d0;delp(n,j,:,:)=0.d0
            enddo
          endif
        endif
      endif
    enddo

  endif

end subroutine collapse
!-----------------------------------------------------------------  

!subroutine collapse(dt,iflag_coll)
!  implicit none
!  real*8,intent(in) :: dt
!  integer,intent(out) :: iflag_coll
!  real*8 rnd,gama_collapse,gama_reset
!  complex*16 su1
!  integer n,i,j
!
!  i=state
!
!  if(icollapse==1) then
!
!    iflag_coll=0
!    do n=1,nquant
!      if(n.ne.state) then
!        gama_reset=sum((force(n,:)-force(i,:))*dble(delr(n,:)-delr(state,:)))/(2*hbar)
!        !! CAUTION !! !! Assumes delr(n,n,:) is in direction of v(:) !!
!        su1=cdabs((V_k(i)-V_k(n))*vdotd(i,n)*sum((delr(n,:)-delr(state,:))*v))/sum(v*v)
!        gama_collapse=gama_reset-2/hbar*cdabs(su1)
!        gama_collapse=gama_collapse*dt
!        gama_reset=-gama_reset*dt
!        call random_number(rnd)
!
!        if(rnd<gama_collapse) then
!          iflag_coll=1
!          cnt_collapse=cnt_collapse+1
!          if(icollapse==1) then
!            !do j=1,nquant
!            !  if(j.ne.n) ci(j)=ci(j)/dsqrt(1-cdabs(ci(n)**2))
!            !enddo
!            !! Erratum: Landry, Subotnik JCP 137, 229901 (2012)
!            ci(i)=ci(i)/cdabs(ci(i))*dsqrt(cdabs(ci(i))**2+cdabs(ci(n))**2)
!            ci(n)=0.d0
!
!          endif
!        endif
!        if(rnd<gama_collapse.or.rnd<gama_reset) then
!          if(icollapse==1) then
!            delr(n,:)=0.d0
!            delp(n,:)=0.d0
!          endif
!        endif
!      endif
!    enddo
!
!  endif
!
!end subroutine collapse
!!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  !! nflag=0: Writes various variables as a function of time
  !! nflag=1: writes minimal useful information at the start and end of trajectory
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 t1,t2
  real*8 phase

  !call cpu_time(t1)

  if(nflag==0) then
    if(iwrite==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(4es17.7,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),temperature,state
        write(11,'(es15.5$)')curr_time*1.d15
        write(12,'(5f15.5)')curr_time*1.d15,cdabs(ci(1:2))**2,datan2(dimag(ci(1:2)),real(ci(1:2)))*180/pi
        write(13,'(5es15.5)')curr_time*1.d15,vdotd(1,2),dasin(W_overlap(1,2))/dtc,hop_prob_net(3-state),state*1.d0
        write(14,'(6f15.5)')curr_time*1.d15,W_overlap(1,1:2),W_overlap(2,1:2),determinant(W_overlap,nquant)
        !write(15,'(6es15.5)')curr_time*1.d15,delr(1,1,1)*1.d10,delr(2,2,1)*1.d10
        !do i=1,nclass
        !  write(11,'(2es15.5$)')x(i)*1.d10,v(i)
        !enddo
        !write(11,*)
      endif
    endif
  endif

  if(nflag==1) then
    if(iwrite==0)then
      write(10,'(5es15.5)')traj_num,energy/wave_to_J,sum(cdabs(ci)**2),temperature
      write(11,*) traj_num
      write(11,'(es15.5$)')curr_time*1.d15
      do i=1,nclass
!        write(11,'(2es15.5$)')x(i)*1.d10,v(i)
      enddo
      write(11,*)
      write(11,*)
    endif
    if(iwrite==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nst_av))/dfloat(nst_av))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      write(13,*);write(13,*)
      write(14,*);write(14,*)
      write(15,*);write(15,*)
    endif
  endif

  !call cpu_time(t2)
  !tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  !! Writes the final useful output
  implicit none
  integer i,j,i1,k
  real*8 nf,pop_el(2)

  nf=dfloat(n_traj)
  cnt_frust=cnt_frust/nf
  cnt_collapse=cnt_collapse/nf

  pop=pop/nf
  rho=rho/nf
  pop_surf=pop_surf/nf
  pop_amp=pop_amp/nf

!do j=1,nquant
!rho1(:)=rho1(:)+rho(j,j,:)
!rho2(:)=rho2(:)+rho(j,j+nquant/2,:)
! enddo

rho1(:)=rho1(:)+rho(1,1,:)
rho2(:)=rho2(:)+rho(2,2,:)
rho3(:)=rho3(:)+rho(3,3,:)
!rho4(:)=rho4(:)+rho(4,4,:)
!rho5(:)=rho5(:)+rho(5,5,:)
!rho6(:)=rho6(:)+rho(6,6,:)
!rho7(:)=rho7(:)+rho(7,7,:)

  do i=1,(nsteps)/nstep_avg
      write(100,'(16f15.7)')(i-1)*nstep_avg*dtc*1.d15,real(rho1(i)),real(rho2(i)),real(rho3(i))!,real(rho4(i)),real(rho5(i)),real(rho6(i)),real(rho7(i))
  enddo

  write(101,*) Vc/wave_to_J,cnt_frust,cnt_collapse

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j

  if(flag==0) then
    !! position dependant variables only
    call tise
    do i=1,nquant
      do j=1,nquant
        sigma(i,j)=ci(i)*dconjg(ci(j))
      enddo
    enddo
  endif

  if(flag==1) then
    KE_en=0.d0
!    do i=1,nclass
!      KE_en=KE_en+0.5*mass(i)*v(i)*v(i)
!    enddo
do i=1,nsite
do j=1,nbath
        KE_en=KE_en+0.5*mass_b*v(j,i)**2
enddo
enddo

    energy=pot_en+KE_en
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! time independent schrodinger equation
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 Hamil(nbasis,nbasis),ens(nbasis),vect(nbasis,nquant)
  real*8 pot_cl,acc_cl(nbath,nsite),acc_qm(nbath,nsite),dpotcl_dx(nbath,nsite)
  real*8 si_adiab_old(nquant,nbasis)
  real*8 t1,t2,part_func

  !call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)
  Hamil_diab=Hamil
  call diag(Hamil,nbasis,ens,vect,nquant)

  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_prev(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

do j=1,nsite
  do i=1,nbath
    delH_dels_ad(state,state,i,j)=sum(si_adiab(:,state)*matmul(delH_dels(:,:,i,j),si_adiab(:,state)))
  enddo
enddo

  call potential_classical(pot_cl,dpotcl_dx)
  acc_qm=-1.d0/mass_b*delH_dels_ad(state,state,:,:)
  acc_cl=-1.d0/mass_b*dpotcl_dx

  pot_en=pot_cl+ens(state)
  V_k=pot_cl+ens(1:nquant)
  acc=acc_cl+acc_qm

  !call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine compute_delH_dels_ad
  implicit none
  integer i,k,kp,i1,j

  force=0.d0
  pot=0.d0
  do k=1,nquant
    do kp=k,nquant
      do i=1,nsite
       do j=1,nbath
        delH_dels_ad(k,kp,j,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,j,i),si_adiab(:,kp)))
      enddo
      enddo
      force(k,kp,:,:)=-delH_dels_ad(k,kp,:,:)
      force(kp,k,:,:)=-delH_dels_ad(k,kp,:,:)
      delH_dels_ad(kp,k,j,i)=delH_dels_ad(k,kp,j,i)
    enddo
    pot(k,k)=V_k(k)
  enddo

  delF=force
  do i1=1,nquant
    delF(i1,i1,:,:)=delF(i1,i1,:,:)-force(state,state,:,:)
  enddo

end subroutine compute_delH_dels_ad
!-----------------------------------------------------------------  

subroutine compute_dij
  implicit none
  integer i,k,kp

!not used
!  do k=1,nquant-1
!    do kp=k+1,nquant
!      do i=1,nclass
!        d_ij(k,kp,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
!      enddo
!      d_ij(k,kp,:)=d_ij(k,kp,:)/(V_k(kp)-V_k(k))
!      d_ij(kp,k,:)=-d_ij(k,kp,:)
!    enddo
!  enddo

end subroutine compute_dij
!-----------------------------------------------------------------  

subroutine compute_dij_2state(x_hop,k,kp,dp)
  implicit none
  integer,intent(in):: k,kp
  real*8,intent(in):: x_hop(nbath,nsite)
  real*8,intent(out):: dp(nbath,nsite)
  real*8 x_sav(nbath,nsite)
  integer i,j

  x_sav=x
  x=x_hop
  call evaluate_variables(0)

  do i=1,nsite
  do j=1,nbath
    dp(j,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,j,i),si_adiab(:,kp)))
  enddo
  enddo
  dp=dp/(V_k(kp)-V_k(k))

  x=x_sav
  call evaluate_variables(0)

end subroutine compute_dij_2state
!-----------------------------------------------------------------  

subroutine compute_vdotd
  !! T matrix computation
  implicit none
  integer i,j,k
  real*8,dimension(nquant,nquant) :: W,ci_W,si_W
  real*8 A,B,C,D,E
  real*8 Wlj,Wlk

  !Method 1
  !call compute_dij
  !vdotd=0.d0
  !do i=1,nclass
  !  vdotd=vdotd+v(i)*d_ij(:,:,i)
  !enddo

  !Method 2
  ! Meek, Levine, JPCL 5, 2351 (2014). Look at Supp info.
!  do j=1,nquant
!    do k=1,nquant
!      W(j,k)=sum(si_adiab_prev(:,j)*si_adiab(:,k))
!      ci_W(j,k)=dacos(W(j,k))
!      si_W(j,k)=dasin(W(j,k))
!    enddo
!  enddo
!
!  vdotd=0.d0
!  do k=1,nquant-1
!    do j=k+1,nquant
!      A=-sinx_x(ci_W(j,j)-si_W(j,k))
!      B=sinx_x(ci_W(j,j)+si_W(j,k))
!      C=sinx_x(ci_W(k,k)-si_W(k,j))
!      D=sinx_x(ci_W(k,k)+si_W(k,j))
!      Wlj=dsqrt(1.d0-W(j,j)**2-W(k,j)**2)
!      if(Wlj==0.d0.or.nquant==2) then
!        E=0.d0
!      else
!        Wlk=(-W(j,k)*W(j,j)-W(k,k)*W(k,j))/Wlj
!        E=2*dasin(Wlj)/(dasin(Wlj)**2-dasin(Wlk)**2)
!        E=E*(Wlj*Wlk*dasin(Wlj)+dasin(Wlk)*(dsqrt((1-Wlj**2)*(1-Wlk**2))-1.d0))
!      endif
!      vdotd(k,j)=0.5/dtc*(ci_W(j,j)*(A+B)+si_W(k,j)*(C+D)+E)
!      vdotd(j,k)=-vdotd(k,j)
!    enddo
!  enddo

  !Method 3
  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab_prev(:,i)*si_adiab(:,j))
    enddo
  enddo

  if(flag_ortho==1)call orthoganalize(W_overlap,nquant)
  call logm(W_overlap,vdotd,nquant)
  vdotd=vdotd/dtc

end subroutine compute_vdotd
!-----------------------------------------------------------------  

subroutine orthoganalize(mat,n)
  integer,intent(in)::n
  real*8,intent(inout)::mat(n,n)
  real*8 S_mat(n,n)

  S_mat=matmul(transpose(mat),mat)
  call inverse_squareroot(S_mat,n)
  mat=matmul(mat,S_mat)

end subroutine orthoganalize
!-----------------------------------------------------------------  
subroutine setup_parameters
  implicit none
  integer i
  real*8 si_diab(nbasis,2),Vb
  real*8 c_0,c_e
  real*8 delw

  mass=1836.d0*au2kg
  mass_b=1836.d0*au2kg

  omg=omg_B
  omg_max=10d0*omg_B
  w_c=omg_max/10d0
  delw=omg_max/real(nclass)

  do i=1,nclass
    !omg(i)=i*delw
    !omg(i)=-w_c*log((i-0.5d0)/float(nclass))
    !ck(i)=dsqrt(2*spectral(omg(i))*mass(i)*omg(i)*delw/pi)
    !ck(i)=omg(i)*dsqrt(2d0*gamma_B*(mass(1)**2)*w_c/(nclass*pi))
  enddo

!omg_photon=(1.549766d0-1.549d0)*ev_2_wave*wave_to_J
!site_en(2)=(1.35d0-1.549d0)*ev_2_wave*wave_to_J
!site_en(1)=(1.549d0-1.549d0)*ev_2_wave*wave_to_J
!site_en(2)=(1.52d0-1.549d0)*ev_2_wave*wave_to_J
!g_coup=dsqrt(2.d0*mass(1)*lambda_B)*omg_B
!coup_intra=-1.d-3*ev_2_wave*wave_to_J
!coup_photon=9.d-3*ev_2_wave*wave_to_J

omg_photon=20000.d0*wave_to_J
site_en=20000.d0*wave_to_J
coup_intra=10*wave_to_J
!coup_photon=1000.d0*wave_to_J

!call SHO

end subroutine setup_parameters

!-----------------------------------------------------------------  

subroutine setup_quantized_vib
  implicit none
  integer i,j,k1,k2
  real*8 H_dvr(n_dvr,n_dvr),ke_dvr(n_dvr,n_dvr),x_dvr(n_dvr),q,delq
  real*8 ens(n_dvr),vect(n_dvr,n_dvr)

  do i=1,n_dvr
    x_dvr(i)=-1.d-10+2.d-10*(i-1)/real(n_dvr-1)
  enddo
  delq=x_dvr(2)-x_dvr(1)

!  en_sho(:,1)=ens(1:nb_vib)

!write(6,*) en_sho/wave_to_J
!stop
  call compute_KE_matrix_dvr(KE_dvr,n_dvr,delq,mass(1))
  H_dvr=KE_dvr
  do i=1,n_dvr
    q=x_dvr(i)
    H_dvr(i,i)=H_dvr(i,i)+0.5*mass(1)*omg_scaled**2*q**2+g_coup*q
  enddo
  call diag(H_dvr,n_dvr,ens,vect,n_dvr)
  si_sho(:,:,1)=vect(:,1:nb_vib)
  en_sho(:,1)=ens(1:nb_vib)


  H_dvr=KE_dvr
  do i=1,n_dvr
    q=x_dvr(i)
    H_dvr(i,i)=H_dvr(i,i)+0.5*mass(1)*omg_scaled**2*q**2-g_coup*q
  enddo
  call diag(H_dvr,n_dvr,ens,vect,n_dvr)
  si_sho(:,:,2)=vect(:,1:nb_vib)
  en_sho(:,2)=ens(1:nb_vib)

  do k1=1,2
    do i=1,nb_vib
      do j=1,nb_vib
        q_exp(i,j,k1)=sum(si_sho(:,i,k1)*x_dvr*si_sho(:,j,k1))
        qsq_exp(i,j,k1)=sum(si_sho(:,i,k1)*x_dvr*x_dvr*si_sho(:,j,k1))
      enddo
    enddo
  enddo

  do k1=1,2
    do k2=1,2
      do i=1,nb_vib
        do j=1,nb_vib
          sho_overlap(i,j,k1,k2)=sum(si_sho(:,i,k1)*si_sho(:,j,k2))
        enddo
      enddo
    enddo
  enddo

  !do i=1,n_dvr
  !  write(20,*) x_dvr(i)*1.d10,si_sho_1(i,1:2)
  !  write(21,*) x_dvr(i)*1.d10,si_sho_2(i,1:2)
  !enddo
  !stop


end subroutine setup_quantized_vib
!-----------------------------------------------------------------  

subroutine setup_H0
  implicit none
  integer i,j,k1,k2
  integer l,m
  real*8 tmp

  Hamil_diab_0=0.d0

  do k1=1,2
    do i=1,nb_vib
      l=(k1-1)*nb_vib+i
      Hamil_diab_0(l,l)=Hamil_site(k1,k1)+en_sho(i,k1)
    enddo
  enddo

  do k1=1,1
    do k2=k1+1,2
      do i=1,nb_vib
        do j=1,nb_vib
          l=(k1-1)*nb_vib+i
          m=(k2-1)*nb_vib+j
          Hamil_diab_0(l,m)=Hamil_diab_0(l,m)+Hamil_site(k1,k2)*sho_overlap(i,j,k1,k2)
          Hamil_diab_0(m,l)=Hamil_diab_0(l,m)
        enddo
      enddo
    enddo
  enddo

  tmp=sum(ck**2/(2*mass(1)*omg**2))

  !do k1=1,3
  !  l=(k1-1)*nb_vib+1
  !  m=(k1)*nb_vib
  !  Hamil_diab_0(l:m,l:m)=Hamil_diab_0(l:m,l:m)+qsq_exp(:,:,k1)*tmp
  !enddo

end subroutine setup_H0
!-----------------------------------------------------------------  

function spectral(w)
  implicit none
  real*8 spectral,w

!  spectral=lambda_B*2.d0*gamma_B*w*omg_B**2/((w**2-omg_B**2)**2+gamma_B**2*w**2)
   spectral = 2.d0*lambda_B*gamma_B*w/(gamma_B**2+w**2)

end function spectral
!-----------------------------------------------------------------  
subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(nquant,nquant),delV_dels(nquant,nquant,nbath,nsite)
  real*8 H1,H2,H12,H3
  real*8,dimension(nclass)::grad_H1,grad_H2,grad_H12
  integer i,j

H_diab=0d0
delV_dels=0.d0

!x=0.d0

do i=1,nsite
        H_diab(i,i)=site_en(i)
        H_diab(i,i+1)=coup_intra
        H_diab(i+1,i)=coup_intra
   do j=1,nbath
        H_diab(i,i)=H_diab(i,i)-ci_bath(j,i)*x(j,i)+ci_bath(j,i)**2/(2*mass_b*omg_bath(j,i)**2)
        delV_dels(i,i,j,i)=-ci_bath(j,i)
   enddo
enddo

H_diab(nbasis,nbasis)=omg_photon
do i=1,nbasis-1
H_diab(i,nbasis)=coup_photon
H_diab(nbasis,i)=coup_photon
enddo

!!! photosynthetic unit hamiltonian
!H_diab(1,1)=0.d0
!H_diab(2,2)=-233.90d0*wave_to_J
!H_diab(3,3)=-233.90d0*wave_to_J
!H_diab(4,4)=-1064.67d0*wave_to_J
!H_diab(nbasis,nbasis)=411.35d0*wave_to_J

!!FMO system
!H_diab(1,1)=200.d0*wave_to_J
!H_diab(2,2)=320.d0*wave_to_J
!H_diab(3,3)=0.d0
!H_diab(4,4)=600.d0*wave_to_J

!H_diab(1,2)=-87.7d0*wave_to_J;  H_diab(2,1)=H_diab(1,2)
!H_diab(1,3)=5.5d0*wave_to_J;    H_diab(3,1)=H_diab(1,3)
!H_diab(2,3)=30.8d0*wave_to_J;   H_diab(3,2)=H_diab(2,3)
!H_diab(1,4)=150.d0*wave_to_J;   H_diab(4,1)=H_diab(1,4)
!H_diab(2,4)=150.d0*wave_to_J;   H_diab(4,2)=H_diab(2,4)
!H_diab(3,4)=150.d0*wave_to_J;   H_diab(4,3)=H_diab(3,4)

!x=0.d0

!H_diab(1,1)=site_en(1)
!H_diab(2,2)=site_en(2)
!H_diab(3,3)=omg_photon
!H_diab(4,4)=(1.36379d0-1.549d0)*eV_2_wave*wave_to_J

!do i=1,nsite
!        do j=1,nbath
!        H_diab(i,i)=H_diab(i,i)-ci_bath(j,i)*x(j,i)+ci_bath(j,i)**2/(2.d0*mass_b*omg_bath(j,i)**2)
!        delV_dels(i,i,j,i)=-ci_bath(j,i)
!        enddo
!enddo

!H_diab(1,2)=coup_intra
!H_diab(2,1)=coup_intra
!H_diab(1,3)=coup_photon
!H_diab(3,1)=coup_photon
!H_diab(2,3)=coup_photon
!H_diab(3,2)=coup_photon
!H_diab(1,4)=coup_photon
!H_diab(4,1)=coup_photon
!H_diab(2,4)=coup_photon
!H_diab(4,2)=coup_photon

!H_diab(1,2)=-24.19d0*wave_to_J;         H_diab(2,1)=H_diab(1,2)
!H_diab(1,3)=-8.0657d0*wave_to_J;        H_diab(3,1)=H_diab(1,3)
!H_diab(2,3)=-266.169d0*wave_to_J;       H_diab(3,2)=H_diab(2,3)
!H_diab(4,5)=275.848d0*wave_to_J;        H_diab(5,4)=H_diab(4,5)
!H_diab(1,5)=72.59d0*wave_to_J;          H_diab(5,1)=H_diab(1,5)
!H_diab(2,5)=72.59d0*wave_to_J;          H_diab(5,2)=H_diab(2,5)
!H_diab(3,5)=72.59d0*wave_to_J;          H_diab(5,3)=H_diab(3,5)

!H_diab(1,2)=coup_intra; H_diab(2,1)=H_diab(1,2)
!H_diab(3,3)=omg_photon
!H_diab(1,3)=coup_photon; H_diab(3,1)=H_diab(1,3)
!H_diab(2,3)=coup_photon; H_diab(3,2)=H_diab(2,3)


!include DSE terms
if(DSE==1)then

do i=1,nbasis-1
do j=1,nbasis-1
H_diab(i,j)=H_diab(i,j)+coup_photon**2/omg_photon
enddo
enddo

H_diab(nbasis,nbasis)=H_diab(nbasis,nbasis)+nsite*coup_photon**2/omg_photon
!H_diab(1,1)=H_diab(1,1)+coup_photon**2/omg_photon
!H_diab(2,2)=H_diab(2,2)+coup_photon**2/omg_photon
!H_diab(3,3)=H_diab(3,3)+2*coup_photon**2/omg_photon
!H_diab(1,2)=H_diab(1,2)+coup_photon**2/omg_photon
!H_diab(2,1)=H_diab(2,1)+coup_photon**2/omg_photon
endif

!write(201,'(10es)')curr_time*1.d15,H_diab(1,1)/wave_to_J,H_diab(2,2)/wave_to_J,H_diab(3,3)/wave_to_J,H_diab(1,2)/wave_to_J,H_diab(1,3)/wave_to_J,H_diab(2,3)/wave_to_J
!do i=1,nquant
!write(6,'(10es)')H_diab(i,:)/wave_to_J
!enddo
!stop

end subroutine compute_potential

!-----------------------------------------------------------------  

subroutine bath_coupling
implicit none

integer::i,j,nbath1
real*8::lambda1,lambda2

!lambda1=0.d0
!lambda2=0.d0
!nsite=1
!del_bath_freq=(bath_freq_max-bath_freq_min)/real(nbath)
!do i=1,nsite
!   do j=1,nbath
!        omg_bath(j,i)=bath_freq_min+j*del_bath_freq
!        ci_bath(j,i)=dsqrt(spectral(omg_bath(j,i))*mass_b*omg_bath(j,i)*del_bath_freq*2/pi)
!        write(40,*)omg_bath(j,i)/(2*pi*clight),spectral(omg_bath(j,i))/wave_to_J
!        lambda1=lambda1+0.5*ci_bath(j,1)**2/(mass_b*omg_bath(j,1)**2)
!   enddo
!enddo

nbath1=nbath
!nsite=1
do i=1,nsite
do j=1,nbath
        omg_bath(j,i)=gamma_B*dtan((real(j)-0.5)*pi/(2*real(nbath1)))
        ci_bath(j,i)=dsqrt(2*mass_b*lambda_B/(real(nbath1)))*omg_bath(j,i)
!        write(42,*)omg_bath(j,i)/(2*pi*clight),spectral(omg_bath(j,i))/wave_to_J
!        lambda2=lambda2+0.5*ci_bath(j,1)**2/(mass_b*omg_bath(j,1)**2)
enddo
enddo
!write(6,*)lambda2/wave_to_J
!stop


!do j=1,nbath
!lambda1=lambda1+0.5*ci_bath(j,1)**2/(mass_b*omg_bath(j,1)**2)
!lambda2=lambda2+0.5*ci_bath(j,2)**2/(mass_b*omg_bath(j,2)**2)
!enddo

end subroutine
!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,acc_cl)
  implicit none
  real*8,intent(out) :: pot_cl,acc_cl(nbath,nsite)
  integer i,j
  real*8 q1,q3

  pot_cl=0.d0
  acc_cl=0.d0

do i=1,nsite
do j=1,nbath
pot_cl=pot_cl+0.5d0*mass_b*omg_bath(j,i)**2*x(j,i)**2!+ci_bath(j,i)**2/(2*mass_b*omg_bath(j,i)**2)
acc_cl(j,i)=mass_b*omg_bath(j,i)**2*x(j,i)
enddo
enddo

end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine check_acceleration
  !! A test subroutine that compares analytical accelerations with numerical
  !accelerations
  implicit none
  integer i,nflag,j
  real*8 delx,en_old,acc_sav(nbath,nsite)
  real*8 q0,rnd

  delx=1.d-17
  state=1

  do i=1,nsite
  do j=1,nbath
    call random_number(rnd)
    x(j,i)=(rnd*2-1.d0)*1.d-10
  enddo
  enddo

  call init_cond

  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(6,*) "delx=",delx
  write(6,*)

  do i=1,nsite
  do j=1,nbath
      x(j,i)=x(j,i)+delx
      call evaluate_variables(0)
      acc(j,i)=-(pot_en-en_old)/delx/mass_b
      write(6,*)"Analytical acceleration =",acc_sav(j,i)
      write(6,*)"Numerical acceleration  =",acc(j,i)
      write(6,*)"Error =",(acc(j,i)-acc_sav(j,i))/acc(j,i)*100.d0
      write(6,*)
      x(j,i)=x(j,i)-delx
  enddo
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

subroutine stochastic_force(delr,delv,dt)
  !! stoachastic forces for langevin equation
  !! Not used for the Holstein model results 
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass

    sig_r=dt*dsqrt(kb*temperature/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

function commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  complex*16 tmp
  integer j,k

  if(iflag==0) commute=matmul(A,B)-matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do j=1,nquant
      do k=1,nquant
        commute(j,k)=B(j,k)*(A(j,j)-A(k,k))
      enddo
    enddo
  endif

  if(iflag==2) then
    !! Assume A is tridiagonal, with a_ii=0, and a_ij=-a_ji (a is assumed to be d_ij)
    do j=1,nquant
      do k=1,nquant
        tmp=0.d0
        if(j<nquant) tmp=tmp+A(j,j+1)*B(j+1,k)
        if(j>1) tmp=tmp-A(j-1,j)*B(j-1,k)
        if(k>1) tmp=tmp-A(k-1,k)*B(j,k-1)
        if(k<nquant) tmp=tmp+A(k,k+1)*B(j,k+1)
      enddo
    enddo
  endif

end function commute
!-----------------------------------------------------------------  

function anti_commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 anti_commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  integer i,j

  if(iflag==0) anti_commute=matmul(A,B)+matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do i=1,nquant
      do j=1,nquant
       anti_commute(i,j)=B(i,j)*(A(i,i)+A(j,j))
      enddo
    enddo
  endif

end function anti_commute
!-----------------------------------------------------------------  

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

subroutine logm(mat,log_mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(in):: mat(n,n)
  real*8,intent(out):: log_mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
  enddo

  log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm
!-----------------------------------------------------------------  

subroutine inverse_squareroot(mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(inout):: mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=1.d0/t(i,i)**0.5d0
  enddo

  mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine inverse_squareroot
!-----------------------------------------------------------------  

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: nold
  complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
  real*8,intent(in) :: mat(n,n)
  complex*16,intent(out) :: T(n,n)
  complex*16,allocatable,intent(inout):: cwork(:)
  real*8 rwork(n)
  complex*16 mat_c(n,n)

  integer lwork
  logical:: select
  logical bwork(n)
  integer sdim,info,AllocateStatus

  T=mat

  info=0
  sdim=0

  if(nold.ne.n .or. .not.allocated(cwork)) then
  !if(nold.ne.n) then
    lwork=-1
    if(allocated(cwork))deallocate(cwork)
    allocate(cwork(n))
    call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
    lwork=int(cwork(1))
    deallocate(cwork)
    allocate(cwork(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(cwork)
  call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine schur
!---------------------------------------------------------- 

REAL FUNCTION determinant(matrix, n)
    !!http://web.hku.hk/~gdli/UsefulFiles/Example-Fortran-program.html
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                determinant= 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    determinant= l
    DO i = 1, n
        determinant= determinant* matrix(i,i)
    END DO
    
END FUNCTION determinant
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!---------------------------------------------------------- 

subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!---------------------------------------------------------- 
subroutine SHO
  !! find eigen functions/ eigen values of H11
  integer :: i,j,k

 !x_grid_max=32d0*g_coup/(dsqrt(2d0)*mass(1)*omg_B**2)
 !x_grid_min=-32d0*g_coup/(dsqrt(2d0)*mass(1)*omg_B**2)
 x_grid_max=1d-10
 x_grid_min=-1d-10

!write(*,*) wave_to_j,0.5d0*hbar*omg_B/wave_to_j,(g_coup**2/(2d0*mass(1)*omg_B**2))/wave_to_j

dx=(x_grid_max-x_grid_min)/real(g-1)

do k=1,g
x_grid(k)=x_grid_min+(k-1)*dx
enddo

do i=1,g
!write(*,*) x_grid(i),hbar*omg_B/wave_to_j
do j=1,g
if (i==j) then
sho_hamil_1(i,j)=0.5d0*mass(1)*(omg_B**2)*(x_grid(i)**2)+(g_coup*x_grid(i))+(((hbar*pi)**2)/(6.d0*mass(1)*(dx**2)))+sum(((ck**2)/(2d0*mass(1)*omg**2)))*(x_grid(i)**2)
else
sho_hamil_1(i,j)=(hbar**2)*((-1)**(i-j))/(mass(1)*(dx**2)*((i-j)**2))
endif
!write(1,*) sho_hamil_1(i,j)
enddo
enddo

call diag(sho_hamil_1,g,ens1,vect1,g)
!write(*,*) ens1/wave_to_j

do i=1,g
x_grid(i)=x_grid_min+(i-1)*dx
do j=1,g
if (i==j) then
sho_hamil_2(i,j)=0.5d0*mass(1)*(omg_B**2)*(x_grid(i)**2)-(g_coup*x_grid(i))+(((hbar*pi)**2)/(6.d0*mass(1)*(dx**2)))+sum(((ck**2)/(2d0*mass(1)*omg**2)))*(x_grid(i)**2)
else
sho_hamil_2(i,j)=(hbar**2)*((-1)**(i-j))/(mass(1)*(dx**2)*((i-j)**2))
endif
enddo
enddo

call diag(sho_hamil_2,g,ens2,vect2,g)
!write(*,*) ens2/wave_to_j

do i=1,g
x_grid(i)=x_grid_min+(i-1)*dx
do j=1,g
if (i==j) then
sho_hamil_0(i,j)=0.5d0*mass(1)*(omg_B**2)*(x_grid(i)**2)+(((hbar*pi)**2)/(6.d0*mass(1)*(dx**2)))+sum(((ck**2)/(2d0*mass(1)*omg**2)))*(x_grid(i)**2)
else
sho_hamil_0(i,j)=(hbar**2)*((-1)**(i-j))/(mass(1)*(dx**2)*((i-j)**2))
endif
enddo
enddo

call diag(sho_hamil_0,g,ens0,vect0,g)
!write(*,*)ens0/wave_to_j
!do i=1,g
!write(30,*)x_grid(i)*1d10,vect0(i,1),vect1(i,1) 
!enddo

!write(*,*) ens0/wave_to_j
!write(*,*) ''
!write(*,*) ens1/wave_to_j
!write(*,*) ''
!write(*,*) ens2/wave_to_j
!stop

do i=1,g
do j=1,g
expec_q11(i,j)=sum(vect1(:,i)*x_grid*vect1(:,j))
expec_q22(i,j)=sum(vect2(:,i)*x_grid*vect2(:,j))
expec_q12(i,j)=sum(vect1(:,i)*x_grid*vect2(:,j))
expec_q21(i,j)=sum(vect2(:,i)*x_grid*vect1(:,j))
expec_sqr_q11(i,j)=sum(vect1(:,i)*(x_grid**2)*vect1(:,j))
expec_sqr_q22(i,j)=sum(vect2(:,i)*(x_grid**2)*vect2(:,j))
expec_sqr_q12(i,j)=sum(vect1(:,i)*(x_grid**2)*vect2(:,j))
expec_sqr_q21(i,j)=sum(vect2(:,i)*(x_grid**2)*vect1(:,j))
expec11(i,j)=sum(vect1(:,i)*vect1(:,j))
expec22(i,j)=sum(vect2(:,i)*vect2(:,j))
expec12(i,j)=sum(vect1(:,i)*vect2(:,j))
expec21(i,j)=sum(vect2(:,i)*vect1(:,j))
enddo
enddo
!stop
!write(9,*)'sho_overlap'
!write(9,*) expec11(1,1:g)
!write(9,*)'q_exp_1'
!write(9,*) expec_q11(1,1:g)
!write(9,*)'q_exp_2'
!write(9,*) expec_q22(1,1:g)
!write(9,*)'qsq_exp_1'
!write(9,*) expec_sqr_q11(1,1:g)
!!write(9,*)'qsq_exp_2'
!write(9,*) expec_sqr_q22(1,1:g)

!stop

!do i=1,g
!do j=1,g
!write(4,*) expec11(i,j)
!enddo
!write(4,*) '' 
!enddo

!stop
!write(1,*) ens1/wave_to_j,ens2/wave_to_j

!do i=1,g
!write(200,*)x_grid(i),vect1(i,1),vect2(i,1) 
!write(201,*)x_grid(i),vect1(i,2),vect2(i,2) 
!write(202,*)x_grid(i),vect1(i,3),vect2(i,3) 
!write(203,*)x_grid(i),vect1(i,4),vect2(i,4) 
!write(204,*)x_grid(i),vect1(i,5),vect2(i,5) 
!enddo

!stop

 !! Construct Hamil_11, and then use diag() to get eigen values/functions.
  !! epsilon 1 = Hamil_site(1,1), epsilon 2 =Hamil_site(2,2), Vc=Hamil_site(1,2)
  !! m = mass(1)
  !! c = ck(1), omega = omg(1)
end subroutine SHO
!----------------------------------------------------
End Module mod_afssh
