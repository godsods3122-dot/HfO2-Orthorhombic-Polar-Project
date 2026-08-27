!------------------------------------------------------------------------------
!> Gonze-Lee long-range dipole-dipole term, ported from phonopy's C
!> implementation (phonopy/c/dynmat.c: dym_get_recip_dipole_dipole and
!> dym_get_recip_dipole_dipole_q0).
!>
!> Rationale: the short-range hr.dat produced from phonopy force constants has
!> phonopy's Gonze long-range term subtracted. The term added back here must be
!> the *same* partition, otherwise the two do not cancel. phonopy differs from
!> the QE rgd_blk lineage in three ways that matter:
!>   1. the structure phase uses G only (C-type), not G+q (D-type),
!>   2. the Ewald screening Lambda is derived from the cell and the G cutoff,
!>   3. the q=0 sum rule term dd_q0 is Hermitian-symmetrised in (alpha,beta).
!>
!> All internal work here is done in Angstrom / eV, matching phonopy, and the
!> result is converted to Simphony's internal units at the end.
!------------------------------------------------------------------------------
subroutine long_range_phonon_interaction_gonze(k, mat2)
   !> Bulk wrapper: k is fractional in the Origin_cell basis.
   use para
   implicit none
   real(dp), intent(in) :: k(3)
   complex(dp), intent(out) :: mat2(Num_wann, Num_wann)
   real(dp) :: cellA(3,3), posA(3,Origin_cell%Num_atoms), Bm(3,3), qA(3), dt
   integer :: i

   cellA(1,:) = Origin_cell%Rua/Angstrom2atomic
   cellA(2,:) = Origin_cell%Rub/Angstrom2atomic
   cellA(3,:) = Origin_cell%Ruc/Angstrom2atomic
   do i=1,Origin_cell%Num_atoms
      posA(:,i) = Origin_cell%Atom_position_cart(:,i)/Angstrom2atomic
   end do
   call inv3(cellA, Bm, dt)
   qA(:) = k(1)*Bm(:,1) + k(2)*Bm(:,2) + k(3)*Bm(:,3)
   call gonze_dd_core(qA, cellA, posA, Origin_cell%Num_atoms, mat2)
end subroutine long_range_phonon_interaction_gonze


!> Cartesian-q wrapper matching the argument style of
!> long_range_phonon_interaction: q and rec_lat are in 2*pi/alat units,
!> tau in alat units, alat_bohr is that cell's cell_parameters(1).
subroutine long_range_phonon_interaction_gonze_cart(q, mat2, tau, rec_lat, natoms, alat_bohr)
   use para
   implicit none
   integer, intent(in) :: natoms
   real(dp), intent(in) :: q(3), tau(3,natoms), rec_lat(3,3), alat_bohr
   complex(dp), intent(out) :: mat2(Num_wann, Num_wann)
   real(dp) :: cellA(3,3), posA(3,natoms), Rm(3,3), dt, alatA, qA(3)
   integer :: i, j

   alatA = alat_bohr/Angstrom2atomic
   !> rec_lat = inv(cell/alat)^T  =>  cell/alat = inv(rec_lat)^T
   call inv3(rec_lat, Rm, dt)
   do i=1,3
      do j=1,3
         cellA(i,j) = Rm(j,i)*alatA
      end do
   end do
   do i=1,natoms
      posA(:,i) = tau(:,i)*alatA
   end do
   qA(:) = q(:)/alatA
   call gonze_dd_core(qA, cellA, posA, natoms, mat2)
end subroutine long_range_phonon_interaction_gonze_cart


!> 3x3 inverse
subroutine inv3(A, Ai, det)
   use para
   implicit none
   real(dp), intent(in) :: A(3,3)
   real(dp), intent(out) :: Ai(3,3), det
   det =  A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
   Ai(1,1)= (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
   Ai(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
   Ai(3,1)= (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
   Ai(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
   Ai(2,2)= (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
   Ai(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
   Ai(1,3)= (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
   Ai(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det
   Ai(3,3)= (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
end subroutine inv3


subroutine gonze_dd_core(qc_in, cellA_in, posA_in, natoms, mat2)
   use para
   implicit none
   integer, intent(in) :: natoms
   real(dp), intent(in) :: qc_in(3), cellA_in(3,3), posA_in(3,natoms)
   complex(dp), intent(out) :: mat2(Num_wann, Num_wann)

   integer :: i, j, a, b, ig, m1, m2, m3, grad, nG, cnt
   real(dp) :: cellA(3,3), Bmat(3,3), det, Vol, Gc, Lam, GeG, trEps
   real(dp) :: Gv(3), qc(3), phase, minnorm, vtest(3)
   real(dp) :: qdir(3), tol
   complex(dp) :: cfac

   real(dp), allocatable, save :: Glist(:,:), posA(:,:)
   complex(dp), allocatable, save :: ddq0(:,:,:)
   real(dp), save :: sLam, sfac, cellcache(3,3)
   integer, save :: snG, natcache = -1
   logical, save :: initialised = .false.

   complex(dp), allocatable :: ddp(:,:,:,:), ddb(:,:,:,:)

   tol = 1d-5
   cellA = cellA_in
   qc = qc_in

   !> Rebuild the cached G list / dd_q0 whenever the cell or atom count changes
   if (initialised) then
      if (natcache /= natoms .or. maxval(abs(cellcache-cellA)) > 1d-10) then
         initialised = .false.
         if (allocated(Glist)) deallocate(Glist)
         if (allocated(posA))  deallocate(posA)
         if (allocated(ddq0))  deallocate(ddq0)
      end if
   end if

   !--------------------------------------------------------------------------
   ! One-off setup: cell in Angstrom, G list, Lambda, dd_q0
   !--------------------------------------------------------------------------
   if (.not. initialised) then
      call inv3(cellA, Bmat, det)
      Vol = abs(det)
      allocate(posA(3,natoms))
      posA = posA_in
      cellcache = cellA
      natcache = natoms

      !> G_cutoff = (3*num_G_points/(4 pi V))^(1/3) with num_G_points = 300
      Gc = (3d0*300d0/(4d0*Pi)/Vol)**(1d0/3d0)

      !> Lambda from exp(-GeG/(4 Lambda^2)) = 1e-10 at the cutoff
      trEps = Diele_Tensor(1,1)+Diele_Tensor(2,2)+Diele_Tensor(3,3)
      GeG = Gc*Gc*trEps/3d0
      Lam = sqrt(-GeG/4d0/log(1d-10))
      sLam = Lam

      !> smallest non-zero reciprocal vector length, to size the search box
      minnorm = 1d30
      do m1=-1,1
        do m2=-1,1
          do m3=-1,1
            if (m1==0 .and. m2==0 .and. m3==0) cycle
            vtest(:) = m1*Bmat(:,1) + m2*Bmat(:,2) + m3*Bmat(:,3)
            if (sqrt(dot_product(vtest,vtest)) < minnorm) minnorm = sqrt(dot_product(vtest,vtest))
          end do
        end do
      end do
      grad = int(Gc/minnorm) + 2

      !> count then store the G vectors inside the cutoff sphere
      cnt = 0
      do m1=-grad,grad
        do m2=-grad,grad
          do m3=-grad,grad
            Gv(:) = m1*Bmat(:,1) + m2*Bmat(:,2) + m3*Bmat(:,3)
            if (dot_product(Gv,Gv) < Gc*Gc) cnt = cnt + 1
          end do
        end do
      end do
      nG = cnt
      snG = nG
      allocate(Glist(3,nG))
      cnt = 0
      do m1=-grad,grad
        do m2=-grad,grad
          do m3=-grad,grad
            Gv(:) = m1*Bmat(:,1) + m2*Bmat(:,2) + m3*Bmat(:,3)
            if (dot_product(Gv,Gv) < Gc*Gc) then
               cnt = cnt + 1
               Glist(:,cnt) = Gv(:)
            end if
          end do
        end do
      end do

      !> factor = nac_factor * 4 pi / V , nac_factor = 14.399652 eV*Angstrom
      sfac = 14.399652d0*4d0*Pi/Vol

      if (.true.) then
         write(*,'(a)')      ' >>> Gonze (phonopy-compatible) LO-TO term enabled'
         write(*,'(a,f12.6)')'     G_cutoff (1/Ang) : ', Gc
         write(*,'(a,i12)')  '     number of G      : ', nG
         write(*,'(a,f12.6)')'     Lambda           : ', Lam
      end if

      !--------------------------------------------------------------------
      ! dd_q0 : the same sum evaluated at q=0, summed over the second atom,
      !         then Hermitian-symmetrised in (alpha,beta)
      !--------------------------------------------------------------------
      allocate(ddp(3,natoms,3,natoms), ddb(3,natoms,3,natoms))
      call gonze_sum(ddp, Glist, nG, natoms, (/0d0,0d0,0d0/), .false., &
                     (/0d0,0d0,0d0/), posA, Lam, tol)
      call gonze_borns(ddb, ddp, natoms)
      allocate(ddq0(3,3,natoms))
      ddq0 = dcmplx(0d0,0d0)
      do i=1,natoms
         do a=1,3
            do b=1,3
               do j=1,natoms
                  ddq0(a,b,i) = ddq0(a,b,i) + ddb(a,i,b,j)
               end do
            end do
         end do
      end do
      do i=1,natoms
         do a=1,3
            do b=a,3
               cfac = (ddq0(a,b,i) + conjg(ddq0(b,a,i)))/2d0
               ddq0(a,b,i) = cfac
               ddq0(b,a,i) = conjg(cfac)
            end do
         end do
      end do
      deallocate(ddp, ddb)
      initialised = .true.
   end if

   !--------------------------------------------------------------------------
   ! dd at the requested q
   !--------------------------------------------------------------------------
   call inv3(cellA, Bmat, det)
   if (LOTO_qdir_run_set) then
      qdir(:) = LOTO_qdir_run(1)*Bmat(:,1) + LOTO_qdir_run(2)*Bmat(:,2) + LOTO_qdir_run(3)*Bmat(:,3)
   else
      qdir(:) = LOTO_qdir(1)*Bmat(:,1) + LOTO_qdir(2)*Bmat(:,2) + LOTO_qdir(3)*Bmat(:,3)
   endif

   allocate(ddp(3,natoms,3,natoms), ddb(3,natoms,3,natoms))
   call gonze_sum(ddp, Glist, snG, natoms, qc, .true., qdir, posA, sLam, tol)
   call gonze_borns(ddb, ddp, natoms)

   do i=1,natoms
      do a=1,3
         do b=1,3
            ddb(a,i,b,i) = ddb(a,i,b,i) - ddq0(a,b,i)
         end do
      end do
   end do

   !> phonopy returns eV/Angstrom^2 (mass unweighted). Convert to the internal
   !> unit of HmnR, which is THz^2 * eV2Hartree^2 (still mass unweighted here,
   !> the caller divides by sqrt(M_i M_j)).
   !> phonopy builds this term with the structure phase carrying G only
   !> (C-type), whereas HmnR from hr.dat is summed with exp(i k.R) alone.
   !> Convert to that gauge before handing it back.
   mat2 = dcmplx(0d0,0d0)
   do i=1,natoms
      do j=1,natoms
         phase = 2d0*Pi*( qc(1)*(posA(1,i)-posA(1,j)) &
                        + qc(2)*(posA(2,i)-posA(2,j)) &
                        + qc(3)*(posA(3,i)-posA(3,j)) )
         cfac = dcmplx(cos(phase), sin(phase))
         do a=1,3
            do b=1,3
               mat2(3*(i-1)+a, 3*(j-1)+b) = ddb(a,i,b,j)*sfac*cfac &
                    *(PhonopytoTHz**2)*(eV2Hartree**2)
            end do
         end do
      end do
   end do

   deallocate(ddp, ddb)

end subroutine gonze_dd_core


!> Structure sum over G of  K_a K_b / (K.eps.K) * exp(-K.eps.K/(4 Lambda^2))
!> times exp(2 pi i (r_i - r_j).G).  Note the phase uses G only, not G+q.
subroutine gonze_sum(ddp, Glist, nG, natoms, qc, use_qdir, qdir, posA, Lam, tol)
   use para
   implicit none
   integer, intent(in) :: nG, natoms
   real(dp), intent(in) :: Glist(3,nG), qc(3), qdir(3), posA(3,natoms), Lam, tol
   logical, intent(in) :: use_qdir
   complex(dp), intent(out) :: ddp(3,natoms,3,natoms)

   integer :: ig, i, j, a, b
   real(dp) :: Gv(3), qK(3), dpart, nrm, phase, L2, KKm(3,3), qd(3), dn
   complex(dp) :: ph(natoms)

   L2 = 4d0*Lam*Lam
   ddp = dcmplx(0d0,0d0)

   do ig=1,nG
      Gv(:) = Glist(:,ig)
      qK(:) = Gv(:) + qc(:)
      nrm = sqrt(dot_product(qK,qK))
      if (nrm < tol) then
         if (.not. use_qdir) cycle
         dn = sqrt(dot_product(qdir,qdir))
         if (dn < 1d-30) cycle
         qd(:) = qdir(:)/dn
         dpart = qd(1)*(Diele_Tensor(1,1)*qd(1)+Diele_Tensor(1,2)*qd(2)+Diele_Tensor(1,3)*qd(3)) &
               + qd(2)*(Diele_Tensor(2,1)*qd(1)+Diele_Tensor(2,2)*qd(2)+Diele_Tensor(2,3)*qd(3)) &
               + qd(3)*(Diele_Tensor(3,1)*qd(1)+Diele_Tensor(3,2)*qd(2)+Diele_Tensor(3,3)*qd(3))
         do a=1,3
            do b=1,3
               KKm(a,b) = qd(a)*qd(b)/dpart
            end do
         end do
      else
         dpart = qK(1)*(Diele_Tensor(1,1)*qK(1)+Diele_Tensor(1,2)*qK(2)+Diele_Tensor(1,3)*qK(3)) &
               + qK(2)*(Diele_Tensor(2,1)*qK(1)+Diele_Tensor(2,2)*qK(2)+Diele_Tensor(2,3)*qK(3)) &
               + qK(3)*(Diele_Tensor(3,1)*qK(1)+Diele_Tensor(3,2)*qK(2)+Diele_Tensor(3,3)*qK(3))
         do a=1,3
            do b=1,3
               KKm(a,b) = qK(a)*qK(b)/dpart*exp(-dpart/L2)
            end do
         end do
      end if

      do i=1,natoms
         phase = 2d0*Pi*dot_product(posA(:,i), Gv(:))
         ph(i) = dcmplx(cos(phase), sin(phase))
      end do

      do i=1,natoms
         do j=1,natoms
            do a=1,3
               do b=1,3
                  ddp(a,i,b,j) = ddp(a,i,b,j) + KKm(a,b)*ph(i)*conjg(ph(j))
               end do
            end do
         end do
      end do
   end do

end subroutine gonze_sum


!> dd(a,i,b,j) = sum_{k,l} Z(i,k,a) * ddp(k,i,l,j) * Z(j,l,b)
subroutine gonze_borns(ddb, ddp, natoms)
   use para
   implicit none
   integer, intent(in) :: natoms
   complex(dp), intent(in) :: ddp(3,natoms,3,natoms)
   complex(dp), intent(out) :: ddb(3,natoms,3,natoms)
   integer :: i, j, a, b, kc, lc
   complex(dp) :: s

   ddb = dcmplx(0d0,0d0)
   do i=1,natoms
      do j=1,natoms
         do a=1,3
            do b=1,3
               s = dcmplx(0d0,0d0)
               do kc=1,3
                  do lc=1,3
                     s = s + Born_Charge(i,kc,a)*ddp(kc,i,lc,j)*Born_Charge(j,lc,b)
                  end do
               end do
               ddb(a,i,b,j) = s
            end do
         end do
      end do
   end do

end subroutine gonze_borns
