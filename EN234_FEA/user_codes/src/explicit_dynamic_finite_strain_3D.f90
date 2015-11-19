!     Subroutines for explicit dynamics with finte strain
!     element identifier == 1010
!==========================SUBROUTINE explicit_dynamic_finite_strain_3D ==============================
subroutine explicit_dynamic_finite_strain_3D(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment
    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points, kint, ii, jj, kk, nsvp ! n_state_variables of each integration point

    real (prec)  ::  dstrain(3,3)             ! Strain vector
    real (prec)  ::  tau(3,3), stress(6)                         ! Kirchoff stress tensor and vector
    real (prec)  ::  Bmat(6,length_dof_array), Bbar(6,length_dof_array), tempmatrix1(6,length_dof_array)
    real (prec)  ::  dNbardy(1:n_nodes,1:3), dNdy(1:n_nodes,1:3)
    real (prec)  ::  dF(3,3), F(3,3), Finv(3,3), delta(3,3), J, eta, deta, el_vol
    real (prec)  ::  dL(3,3), dLbar(3,3), dW(3,3), dRot(3,3)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  dof_tot(length_dof_array/3, 3), dof_inc(length_dof_array/3, 3)
    real (prec)  ::  tempL, tempJ, tempInv(3,3)
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    element_deleted = .true.
    x = reshape(element_coords,(/3,length_coord_array/3/))
    dof_tot = transpose( reshape( dof_total + 0.5D0*dof_increment, (/3,length_coord_array/3/) ) )
    dof_inc = transpose( reshape( dof_increment, (/3,length_coord_array/3/) ) )
    delta = reshape((/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), shape(delta))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    el_vol = 0.D0
    eta = 0.D0
    deta = 0.D0
    dNbardy = 0.D0
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! find F, dF
        dF = 0.D0
        F = 0.D0
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, n_nodes
                    dF(ii,jj) = dF(ii,jj) + dof_inc(kk,ii)*dNdx(kk,jj)
                    F(ii,jj) = F(ii,jj) + dof_tot(kk,ii)*dNdx(kk,jj)
                end do
                F(ii,jj) = F(ii,jj) + delta(ii,jj)
            end do
        enddo
        ! find Finv, J
        call invert_small(F,Finv,J)
        ! find dNdy
        dNdy = 0.D0
        do ii = 1, 3
            do jj = 1, 3
                dNdy(1:n_nodes,ii) = dNdy(1:n_nodes,ii) + dNdx(1:n_nodes,jj)*Finv(jj,ii)
            end do
        end do
        ! find dNbardy
        dNbardy(1:n_nodes,1:3) = dNbardy(1:n_nodes,1:3) + J*dNdy(1:n_nodes,1:3)*w(kint)*determinant
        ! find el_vol
        el_vol = el_vol + w(kint)*determinant
        ! find dL
        dL = matmul(dF,Finv)
        ! find eta, deta
        eta = eta + J*w(kint)*determinant
        deta = deta + J*(dL(1,1)+dL(2,2)+dL(3,3))*w(kint)*determinant
    end do
    eta = eta/el_vol
    deta = deta/eta/el_vol
    dNbardy = dNbardy/eta/el_vol

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! find dF, F
        dF = 0.D0
        F = 0.D0
        do ii = 1, 3
            do jj = 1, 3
                do kk = 1, n_nodes
                    dF(ii,jj) = dF(ii,jj) + dof_inc(kk,ii)*dNdx(kk,jj)
                    F(ii,jj) = F(ii,jj) + dof_tot(kk,ii)*dNdx(kk,jj)
                end do
                F(ii,jj) = F(ii,jj) + delta(ii,jj)
            end do
        enddo
        call invert_small(F,Finv,J)
        ! find dNdy
        dNdy = 0.D0
        do ii = 1, 3
            do jj = 1, 3
                dNdy(1:n_nodes,ii) = dNdy(1:n_nodes,ii) + dNdx(1:n_nodes,jj)*Finv(jj,ii)
            end do
        end do
        ! find dLbar
        tempL = 0.D0
        do ii = 1,3
            do jj = 1,3
                tempL = tempL + dF(ii,jj)*Finv(jj,ii)
            end do
        end do
        dLbar = matmul(dF,Finv) + (deta-tempL)/3.D0*delta
        ! find dstrain
        dstrain = (dLbar+transpose(dLbar))/2.D0
        ! find dW
        dW = (dLbar-transpose(dLbar))/2.D0
        ! find dR
        call invert_small(delta-0.5D0*dW,tempInv,tempJ)
        dRot = matmul(tempInv,delta+0.5D0*dW)
	    ! find tau! total 8 state variables of each integration point
        nsvp = 8
        call gurson(element_properties,n_properties,nsvp, &
                    initial_state_variables(1+(kint-1)*nsvp:kint*nsvp),dstrain,dRot, &
                    updated_state_variables(1+(kint-1)*nsvp:kint*nsvp),stress)
        if (updated_state_variables(nsvp*kint) .lt. element_properties(13)) then
            element_deleted = .false.
        end if
        tau = reshape((/ stress(1), stress(4), stress(5), &
                         stress(4), stress(2), stress(6), &
                         stress(5), stress(6), stress(3) /),shape(tau))
        ! find Bbar
        Bmat = 0.D0
        Bmat(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bmat(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bmat(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        Bmat(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        Bmat(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        Bmat(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        Bmat(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        Bmat(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        Bmat(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        tempmatrix1 = 0.D0
        Bbar = 0.D0
        do ii = 1, n_nodes
            tempmatrix1(1:3,3*ii-2) = dNbardy(ii,1) - dNdy(ii,1)
            tempmatrix1(1:3,3*ii-1) = dNbardy(ii,2) - dNdy(ii,2)
            tempmatrix1(1:3,3*ii-0) = dNbardy(ii,3) - dNdy(ii,3)
        end do
        Bbar = Bmat + 1.D0/3.D0*tempmatrix1
	    ! find element force vector
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant

    end do
    return
end subroutine explicit_dynamic_finite_strain_3D

!==========================SUBROUTINE gurson ==============================
subroutine gurson(element_properties,n_properties,n_state_variables,initial_state_variables,dstrain,dRot, &
    updated_state_variables,stress)
    use Types
    use ParamIO
    use Globals, only : TIME, DTIME
    use Element_Utilities, only : rotatesymvec

    implicit none

    integer, intent( in ) :: n_properties
    integer, intent( in ) :: n_state_variables
    real (prec), intent( in ) :: element_properties(n_properties)
    real (prec), intent( in ) :: initial_state_variables(n_state_variables)
    real (prec), intent( in ) :: dstrain(3,3)
    real (prec), intent( in ) :: dRot(3,3)
    real (prec), intent( out ) :: stress(6)
    real (prec), intent( out ) :: updated_state_variables(n_state_variables)

    real (prec) :: taun(3,3), taunp1(3,3), tauDn(3,3), Sstar(3,3), de(3,3), delta(3,3)
    real (prec) :: Vfn, Vfnp1, ematrixn, ematrixnp1, dematrix, stress0(6), stress1(6)
    real (prec) :: tempmat1(3,3)
    real (prec) :: phi2, se, sestar, p, pstar, fstar, fFbar
    real (prec) :: dee, depsnu, dfidse, dfidp
    real (prec) :: E, xnu, Y, epsdot0, m, q1, q2, q3, fN, epsN, sN, fc, fF              ! Material properties
    real (prec) :: a, b, c, h, g, deps, tol
    integer :: ii, jj

    tol = 1.D-16
    delta = reshape((/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0 /), shape(delta))

    E = element_properties(1)
    xnu = element_properties(2)
    Y = element_properties(3)
    epsdot0 = element_properties(4)
    m = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fN = element_properties(9)
    sN = element_properties(10)
    epsN = element_properties(11)
    fc = element_properties(12)
    fF = element_properties(13)

    ! extract n-step state variables
    stress = initial_state_variables(1:6)
    ematrixn = initial_state_variables(7)
    Vfn = initial_state_variables(8)
    taun = reshape((/ stress(1), stress(4), stress(5), &
                      stress(4), stress(2), stress(6), &
                      stress(5), stress(6), stress(3) /),shape(taun))
    h = 3.D0*E/2/(1+xnu)
    g = E/3.D0/(1-2*xnu)
    ! find tauDn
    tauDn = taun-(taun(1,1)+taun(2,2)+taun(3,3))/3.D0*delta
    ! find de
    de = dstrain-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.D0*delta

    ! find pstar
    pstar = (taun(1,1)+taun(2,2)+taun(3,3))/3.D0 + g*(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))
    ! find Sstar
    stress0 = (/tauDn(1,1),tauDn(2,2),tauDn(3,3),tauDn(1,2),tauDn(1,3),tauDn(2,3)/)
    stress1 = rotatesymvec(stress0,dRot)
    tempmat1 = reshape((/ stress1(1), stress1(4), stress1(5), &
                          stress1(4), stress1(2), stress1(6), &
                          stress1(5), stress1(6), stress1(3) /),shape(tempmat1))
    Sstar = h/1.5D0*de + tempmat1
    ! find fFbar
    fFbar = (q1 + dsqrt(q1**2-q3))/q3
    ! find fstar
    if (Vfn .lt. fc) then
        fstar = Vfn
    elseif ((fc .le. Vfn) .and. (Vfn .le. fF)) then
        fstar = fc + (fFbar-fc)*(Vfn-fc)/(fF-fc)
    elseif (Vfn .gt. fF) then
        fstar = fFbar
    endif
    a = 2.D0*q1*fstar
    b = 3.D0*q2/2/Y
    c = 1.D0+q3*fstar**2
    deps = epsdot0*DTIME
    ! find sestar
    sestar = 0.D0
    do ii = 1,3
        do jj = 1,3
            sestar = sestar + 3.D0/2.D0*Sstar(ii,jj)**2
        end do
    end do
    sestar = dsqrt(sestar)
    ! find phi
    phi2 = sestar**2/Y**2+a*cosh(b*pstar)-c
    ! find taunp1, Vfnp1, ematrixnp1
    if (phi2 .lt. tol) then
        ! elastci step
!        write(*,'(A20)') 'Elastic step...'
        ! find taunp1 to update tau
        taunp1 = Sstar + pstar*delta
        ! update Vf, ematrix
        ematrixnp1 = ematrixn
        Vfnp1 = Vfn
    else
        ! plastic step
!        write(*,'(A20)') 'Plastic step...'
        ! find dee & depsnu by N-R loop
        call NRloop(dee, depsnu, a, b, c, h, g, sestar, pstar, Y, m, deps)
        ! find se, p
        se = sestar - h*dee
        p = pstar - g*depsnu
        ! find dfidse, dfidp
        dfidse = se/Y**2/dsqrt(se**2/Y**2+a*cosh(b*p)-c)
        dfidp = a*b*sinh(b*p)/2/dsqrt(se**2/Y**2+a*cosh(b*p)-c)
        ! find taunp1 to update tau
        taunp1 = Sstar - h*dee/sestar*Sstar + (pstar-g*depsnu)*delta
        ! update Vf, ematrix
        dematrix = deps*(dsqrt(phi2)**m)*(dfidse*se+dfidp*p/3.D0)/(1.D0-Vfn)/dsqrt((dfidse**2+2.D0/9*dfidp**2))
        ematrixnp1 = ematrixn + dematrix
        Vfnp1 = 1.D0 + (Vfn-1.D0)*exp(-depsnu) + fN*dematrix/sN/dsqrt(2*PI_D)*exp(-0.5D0*((ematrixn-epsN)/sN)**2)
    end if
    ! upstate stress
    stress = (/taunp1(1,1),taunp1(2,2),taunp1(3,3),taunp1(1,2),taunp1(1,3),taunp1(2,3)/)
    ! upstate n+1-step state variables
    updated_state_variables(1:6) = stress
    updated_state_variables(7) = ematrixnp1
    updated_state_variables(8) = Vfnp1
end subroutine gurson

!==========================SUBROUTINE NRloop ==============================
subroutine NRloop(dee, depsnu, a, b, c, h, g, sestar, pstar, Y, m, deps)
    use Types
    use ParamIO
    use Element_Utilities, only : invert_small
    implicit none
    real (prec)  ::  dee, depsnu, a, b, c, h, g, sestar, pstar, Y, m, deps
    real (prec)  ::  se, dsedee, p, dpdepsnu, phi, dfidse, dfidp, dfi2dse2, dfi2dp2, dfi2dsedp
    real (prec)  ::  f1, f2, Df1Ddee, Df1Ddepsnu, Df2Ddee, Df2Ddepsnu
    real (prec)  ::  tempDf(2,2), tempDfinv(2,2), tempf(2), tempfdet, tol, err
    real (prec)  ::  xc(2), xn(2), dx(2), xnorm, dxnorm
    integer  :: iternum

    tol = 1.D-8
    err = 1.D0
    xc = (/0.D0, 0.D0/)
    iternum = 0
    ! NR loop
    do while (err .gt. tol)
        dee = xc(1)
        depsnu = xc(2)
        se = sestar-h*dee
        dsedee = -h
        p = pstar-g*depsnu
        dpdepsnu = -g
        phi = dsqrt(se**2/Y**2+a*cosh(b*p)-c)
        dfidse = se/Y**2/phi
        dfidp = a*b*sinh(b*p)/2.D0/phi
        dfi2dse2 = 1.D0/Y**2/phi - se/Y**2/phi**2*dfidse
        dfi2dp2 = a*b**2*cosh(b*p)/2.D0/phi-a*b*sinh(b*p)/2.D0/phi**2*dfidp
        dfi2dsedp = -se/Y**2/phi**2-dfidp
        f1 = dsqrt(dfidse**2+2.D0/9*dfidp**2)*dee/deps - dfidse*phi**m
        f2 = dsqrt(dfidse**2+2.D0/9*dfidp**2)*depsnu/deps - dfidp*phi**m
        Df1Ddee = 0.5D0/dsqrt(dfidse**2+2.D0/9*dfidp**2)*(2*dfidse*dfi2dse2+4.D0/9*dfidp*dfi2dsedp)*dsedee*dee/deps &
                     + dsqrt(dfidse**2+2.D0/9*dfidp**2)/deps - dfi2dse2*dsedee*phi**m - dfidse**2*m*phi**(m-1)*dsedee
        Df1Ddepsnu = 0.5D0/dsqrt(dfidse**2+2.D0/9*dfidp**2)*(2*dfidse*dfi2dsedp+4.D0/9*dfidp*dfi2dp2)*dpdepsnu*dee/deps &
                        - dfi2dsedp*dpdepsnu*phi**m - dfidse*dfidp*m*phi**(m-1)*dpdepsnu
        Df2Ddee = 0.5D0/dsqrt(dfidse**2+2.D0/9*dfidp**2)*(2*dfidse*dfi2dse2+4.D0/9*dfidp*dfi2dsedp)*dsedee*depsnu/deps &
                     - dfi2dsedp*dsedee*phi**m - dfidse*dfidp*m*phi**(m-1)*dsedee
        Df2Ddepsnu = 0.5D0/dsqrt(dfidse**2+2.D0/9*dfidp**2)*(2*dfidse*dfi2dsedp+4.D0/9*dfidp*dfi2dp2)*dpdepsnu*depsnu/deps &
                        + dsqrt(dfidse**2+2.D0/9*dfidp**2)/deps - dfi2dp2*dpdepsnu*phi**m - dfidp**2*m*phi**(m-1)*dpdepsnu
        tempDf = reshape((/Df1Ddee, Df2Ddee, Df1Ddepsnu, Df2Ddepsnu/),(/2,2/))
        tempf = -(/f1,f2/)
        call invert_small(tempDf,tempDfinv,tempfdet)
        dx = matmul(tempDfinv,tempf)
        xn = xc + dx
        xnorm = dsqrt(xn(1)**2+xn(2)**2)
        dxnorm = dsqrt(dx(1)**2+dx(2)**2)
        err = dxnorm/xnorm
        xc = xn
        iternum = iternum + 1

!        write(*,'(A12, I2, A12, E12.5, A12, E12.5)') 'Iteration #:', iternum, 'Error:', err, 'Tolerance:', tol

    end do
    dee = xc(1)
    depsnu = xc(2)
end subroutine NRloop


!==========================SUBROUTINE fieldvars_explicit_dynamic_finite_strain_3D ==============================
subroutine fieldvars_explicit_dynamic_finite_strain_3D(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,nsvp
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  stress(6)                         ! Kirchoff stress tensor and true stress vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  sdev(6), p, smises, Vf

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        nsvp = 8
        ! Kirchkoff stress
        stress = updated_state_variables(1+(kint-1)*nsvp:kint*nsvp-2)
        Vf = updated_state_variables(nsvp*kint)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1, n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'Vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vf*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
    end do

    return
end subroutine fieldvars_explicit_dynamic_finite_strain_3D

