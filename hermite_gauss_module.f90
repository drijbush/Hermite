Module hgauss_module

  ! hermite gaussian module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Private :: coeff_array
     Real( wp ), Dimension( : ), Allocatable, Private :: t_data
  End Type coeff_array

  Type, Public :: hgauss_coeffs
     Private
     Integer                                            , Private :: L_alloc_max
     Integer                                            , Private :: L1_max
     Integer                                            , Private :: L2_max
     Real( wp )                                         , Private :: a1, a2
     Real( wp )                                         , Private :: p , q
     Real( wp )                                         , Private :: x12
     Type( coeff_array ), Dimension( :, : ), Allocatable, Private :: coeffs
   Contains
     Procedure, Public :: allocate     => hgauss_alloc
     Procedure, Public :: calc_product => hgauss_calc_product
     Procedure, Public :: calc_coeffs  => hgauss_calc_coeffs
  End Type hgauss_coeffs

  Private

  Logical, Parameter, Private :: do_debug = .True.
  
Contains

  Pure Subroutine hgauss_alloc( Eij, L_max, do_init_with_NaNs )
    
    Class( hgauss_coeffs ), Intent( InOut )           :: Eij
    Integer               , Intent( In    )           :: L_max
    Logical               , Intent( In    ), Optional :: do_init_with_NaNs

    Integer :: l1, l2

    ! Do it this way to make sure all members of Eij are deallocated
    ! before reallocation as can not have a polymorphic dummy argument Intent( Out )
    ! if the passed var
    Call with_dealloc( Eij, L_max )

    If( Present( do_init_with_NaNs ) .And. do_debug ) Then
       If( do_init_with_NaNs ) Then
          Do l2 = 0, L_max
             Do l1 = 0, L_max
                Call init_with_NaNs( Eij%coeffs( l1, l2 )%t_data )
             End Do
          End Do
       End If
    End If

  Contains

    Pure Subroutine with_dealloc( Eij, L_max )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Type( hgauss_coeffs ), Intent(   Out ) :: Eij
      Integer              , Intent( In    ) :: L_max

      Integer :: L
      Integer :: l1, l2

      Eij%L_alloc_max = L_max

      Allocate( Eij%coeffs( 0:L_max, 0:L_max ) )

      Do l2 = 0, L_max
         Do l1 = 0, L_max
            L = l1 + l2
            Allocate( Eij%coeffs( l1, l2 )%t_data( 0:L ) )
         End Do
      End Do
      
    End Subroutine with_dealloc

  End Subroutine hgauss_alloc

  Pure Subroutine hgauss_calc_coeffs( Eij, l1, l2, a1, a2, x12, do_init_with_NaNs )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Implicit None

    Class( hgauss_coeffs ), Intent( InOut )           :: Eij
    Integer               , Intent( In    )           :: l1
    Integer               , Intent( In    )           :: l2
    Real( wp )            , Intent( In    )           :: a1
    Real( wp )            , Intent( In    )           :: a2
    Real( wp )            , Intent( In    )           :: x12
    Logical               , Intent( In    ), Optional :: do_init_with_NaNs

    Real( wp ) :: p, q
    Real( wp ) :: xq
    Real( wp ) :: Kab
   
    Integer :: L_max
    Integer :: li, lj
    
    L_max = l1 + l2

    p = a1 + a2
    q = a1 * a2 / ( a1 + a2 )

    xq = x12
    Kab = Exp( - q * xq * xq )

    Eij%L1_max = l1
    Eij%L2_max = l2

    Eij%a1 = a1
    Eij%a2 = a2

    Eij%p = p
    Eij%q = q

    Eij%x12 = xq
    
    Eij%coeffs( 0, 0 )%t_data( 0 ) = Kab ! L1=0, L2=0 case

    Angular_momentum: Select Case( L_max )
    Case( 0 )
       ! Already dealt with
    Case( 1 )
       Select Case( l1 )
       Case( 0 )
          Call hgauss_calc_coeffs_0_1( a1, a2, xq, Eij )
       Case( 1 )
          Call hgauss_calc_coeffs_1_0( a1, a2, xq, Eij )
       End Select
    Case( 2 )
       Select Case( l1 )
       Case( 0 )
          Call hgauss_calc_coeffs_0_2( a1, a2, xq, Eij )
       Case( 1 )
          Call hgauss_calc_coeffs_1_1( a1, a2, xq, Eij )
       Case( 2 )
          Call hgauss_calc_coeffs_2_0( a1, a2, xq, Eij )
       End Select
    Case Default
       Call hgauss_calc_coeffs_default( l1, l2, a1, a2, xq, Eij )
    End Select Angular_momentum

    If( Present( do_init_with_NaNs ) .And. do_debug ) Then
       If( do_init_with_NaNs ) Then
          Do lj = l2 + 1, Ubound( Eij%coeffs, Dim = 2 )
             Do li = l1 + 1, Ubound( Eij%coeffs, Dim = 1 )
                Call init_with_NaNs( Eij%coeffs( li, lj )%t_data )
             End Do
          End Do
       End If
    End If

  Contains

    Pure Subroutine hgauss_calc_coeffs_0_1( a1, a2, xq, Eij )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Implicit None

      Real( wp )           , Intent( In    ) :: a1
      Real( wp )           , Intent( In    ) :: a2
      Real( wp )           , Intent( In    ) :: xq
      Type( hgauss_coeffs ), Intent( InOut ) :: Eij

      Real( wp ) :: p, q
      Real( wp ) :: pfac, qxq
      Real( wp ) :: qxq_over_a2

      p = a1 + a2
      q = a1 * a2 / ( a1 + a2 )

      pfac = 0.5_wp / p

      qxq  = q * xq
      qxq_over_a2 = qxq / a2

      Eij%coeffs( 0, 1 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 0, 0 )%t_data( 0 )
      Eij%coeffs( 0, 1 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )

    End Subroutine hgauss_calc_coeffs_0_1

    Pure Subroutine hgauss_calc_coeffs_1_0( a1, a2, xq, Eij )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Implicit None

      Real( wp )           , Intent( In    ) :: a1
      Real( wp )           , Intent( In    ) :: a2
      Real( wp )           , Intent( In    ) :: xq
      Type( hgauss_coeffs ), Intent( InOut ) :: Eij

      Real( wp ) :: p, q
      Real( wp ) :: pfac, qxq
      Real( wp ) :: qxq_over_a1

      p = a1 + a2
      q = a1 * a2 / ( a1 + a2 )

      pfac = 0.5_wp / p

      qxq  = q * xq
      qxq_over_a1 = qxq / a1

      Eij%coeffs( 1, 0 )%t_data( 0 ) = - qxq_over_a1 * Eij%coeffs( 0, 0 )%t_data( 0 )
      Eij%coeffs( 1, 0 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )

    End Subroutine hgauss_calc_coeffs_1_0

    Pure Subroutine hgauss_calc_coeffs_0_2( a1, a2, xq, Eij )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Implicit None

      Real( wp )           , Intent( In    ) :: a1
      Real( wp )           , Intent( In    ) :: a2
      Real( wp )           , Intent( In    ) :: xq
      Type( hgauss_coeffs ), Intent( InOut ) :: Eij

      Real( wp ) :: p, q
      Real( wp ) :: pfac, qxq
      Real( wp ) :: qxq_over_a2

      p = a1 + a2
      q = a1 * a2 / ( a1 + a2 )

      pfac = 0.5_wp / p

      qxq  = q * xq
      qxq_over_a2 = qxq / a2

      Eij%coeffs( 0, 1 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 0, 0 )%t_data( 0 )
      Eij%coeffs( 0, 1 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )
      
      Eij%coeffs( 0, 2 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 0, 1 )%t_data( 0 ) + &
                                                       Eij%coeffs( 0, 1 )%t_data( 1 ) 
      Eij%coeffs( 0, 2 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 1 )%t_data( 0 ) + &
                                         qxq_over_a2 * Eij%coeffs( 0, 1 )%t_data( 1 )
      Eij%coeffs( 0, 2 )%t_data( 2 ) =          pfac * Eij%coeffs( 0, 1 )%t_data( 1 )
      
    End Subroutine hgauss_calc_coeffs_0_2

    Pure Subroutine hgauss_calc_coeffs_1_1( a1, a2, xq, Eij )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Implicit None

      Real( wp )           , Intent( In    ) :: a1
      Real( wp )           , Intent( In    ) :: a2
      Real( wp )           , Intent( In    ) :: xq
      Type( hgauss_coeffs ), Intent( InOut ) :: Eij

      Real( wp ) :: p, q
      Real( wp ) :: pfac, qxq
      Real( wp ) :: qxq_over_a1
      Real( wp ) :: qxq_over_a2

      p = a1 + a2
      q = a1 * a2 / ( a1 + a2 )

      pfac = 0.5_wp / p

      qxq  = q * xq
      qxq_over_a1 = qxq / a1
      qxq_over_a2 = qxq / a2

      Eij%coeffs( 1, 0 )%t_data( 0 ) = - qxq_over_a1 * Eij%coeffs( 0, 0 )%t_data( 0 )
      Eij%coeffs( 1, 0 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )

      Eij%coeffs( 0, 1 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 0, 0 )%t_data( 0 )
      Eij%coeffs( 0, 1 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )

      Eij%coeffs( 1, 1 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 1, 0 )%t_data( 0 ) + &
                                                       Eij%coeffs( 1, 0 )%t_data( 1 )
      Eij%coeffs( 1, 1 )%t_data( 1 ) =          pfac * Eij%coeffs( 1, 0 )%t_data( 0 ) + &
                                         qxq_over_a2 * Eij%coeffs( 1, 0 )%t_data( 1 )
      Eij%coeffs( 1, 1 )%t_data( 2 ) =          pfac * Eij%coeffs( 1, 0 )%t_data( 1 )

    End Subroutine hgauss_calc_coeffs_1_1

    Pure Subroutine hgauss_calc_coeffs_2_0( a1, a2, xq, Eij )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Implicit None

      Real( wp )           , Intent( In    ) :: a1
      Real( wp )           , Intent( In    ) :: a2
      Real( wp )           , Intent( In    ) :: xq
      Type( hgauss_coeffs ), Intent( InOut ) :: Eij

      Real( wp ) :: p, q
      Real( wp ) :: pfac, qxq
      Real( wp ) :: qxq_over_a1

      p = a1 + a2
      q = a1 * a2 / ( a1 + a2 )

      pfac = 0.5_wp / p

      qxq  = q * xq
      qxq_over_a1 = qxq / a1

      Eij%coeffs( 1, 0 )%t_data( 0 ) = - qxq_over_a1 * Eij%coeffs( 0, 0 )%t_data( 0 )
      Eij%coeffs( 1, 0 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )
      
      Eij%coeffs( 2, 0 )%t_data( 0 ) = - qxq_over_a1 * Eij%coeffs( 1, 0 )%t_data( 0 ) + &
                                                       Eij%coeffs( 1, 0 )%t_data( 1 ) 
      Eij%coeffs( 2, 0 )%t_data( 1 ) =          pfac * Eij%coeffs( 1, 0 )%t_data( 0 ) - &
                                         qxq_over_a1 * Eij%coeffs( 1, 0 )%t_data( 1 )
      Eij%coeffs( 2, 0 )%t_data( 2 ) =          pfac * Eij%coeffs( 1, 0 )%t_data( 1 )
      
    End Subroutine hgauss_calc_coeffs_2_0

    Pure Subroutine hgauss_calc_coeffs_default( l1, l2, a1, a2, xq, Eij )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

      Implicit None

      Integer              , Intent( In    ) :: l1
      Integer              , Intent( In    ) :: l2
      Real( wp )           , Intent( In    ) :: a1
      Real( wp )           , Intent( In    ) :: a2
      Real( wp )           , Intent( In    ) :: xq
      Type( hgauss_coeffs ), Intent( InOut ) :: Eij

      Real( wp ) :: p, q
      Real( wp ) :: pfac, qxq
      Real( wp ) :: qxq_over_a1, qxq_over_a2
      Real( wp ) :: t_fac

      Integer :: L_max, L
      Integer :: i, j, t
      Integer :: ip1, jp1

      L_max = l1 + l2

      p = a1 + a2
      q = a1 * a2 / ( a1 + a2 )

      pfac = 0.5_wp / p

      qxq  = q * xq
      qxq_over_a1 = qxq / a1
      qxq_over_a2 = qxq / a2

      ! Fill in the j = 0 values for all values of i
      ! First the unusual low vaues of l1 where some terms are zero since they just don't exist
      If( l1 >= 1 ) Then
         Eij%coeffs( 1, 0 )%t_data( 0 ) = - qxq_over_a1 * Eij%coeffs( 0, 0 )%t_data( 0 )
         Eij%coeffs( 1, 0 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )
      End If

      If( l1 >= 2 ) Then
         Eij%coeffs( 2, 0 )%t_data( 0 ) = - qxq_over_a1 * Eij%coeffs( 1, 0 )%t_data( 0 ) + &
                                                          Eij%coeffs( 1, 0 )%t_data( 1 ) 
         Eij%coeffs( 2, 0 )%t_data( 1 ) =          pfac * Eij%coeffs( 1, 0 )%t_data( 0 ) - &
                                            qxq_over_a1 * Eij%coeffs( 1, 0 )%t_data( 1 )
         Eij%coeffs( 2, 0 )%t_data( 2 ) =          pfac * Eij%coeffs( 1, 0 )%t_data( 1 )
      End If

      ! Now the more general terms
      Do L = 3, l1

         ip1 = L
         i = ip1 - 1

         t_fac = 1.0_wp

         t = 0

         Eij%coeffs( ip1, 0 )%t_data( t ) = - qxq_over_a1 * Eij%coeffs( i, 0 )%t_data( t     ) + &
                                                    t_fac * Eij%coeffs( i, 0 )%t_data( t + 1 )

         Do t = 1, L - 2
            t_fac = t_fac + 1.0_wp
            Eij%coeffs( ip1, 0 )%t_data( t ) =        pfac * Eij%coeffs( i, 0 )%t_data( t - 1 ) - &
                                               qxq_over_a1 * Eij%coeffs( i, 0 )%t_data( t     ) + &
                                                     t_fac * Eij%coeffs( i, 0 )%t_data( t + 1 )
         End Do

         t = L - 1
         Eij%coeffs( ip1, 0 )%t_data( t ) =         pfac * Eij%coeffs( i, 0 )%t_data( t - 1 ) - &
                                             qxq_over_a1 * Eij%coeffs( i, 0 )%t_data( t     ) 

         t = L
         Eij%coeffs( ip1, 0 )%t_data( t ) =         pfac * Eij%coeffs( i, 0 )%t_data( t - 1 )

      End Do

      ! OK All j = 0 values filled in. Now can fill in the rest
      If( l2 >= 1 ) Then
         Eij%coeffs( 0, 1 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 0, 0 )%t_data( 0 )
         Eij%coeffs( 0, 1 )%t_data( 1 ) =          pfac * Eij%coeffs( 0, 0 )%t_data( 0 )
      End If

      If( L_max >= 2 ) Then
         L = 2
         Do jp1 = Max( 1, L - 1 - l1 + 1 ), Min( L, l2 )
            j = jp1 - 1
            i = L - 1 - j
            Eij%coeffs( i, jp1 )%t_data( 0 ) = + qxq_over_a2 * Eij%coeffs( i, j )%t_data( 0 ) + &
                                                               Eij%coeffs( i, j )%t_data( 1 )
            Eij%coeffs( i, jp1 )%t_data( 1 ) =          pfac * Eij%coeffs( i, j )%t_data( 0 ) + &
                                                 qxq_over_a2 * Eij%coeffs( i, j )%t_data( 1 )
            Eij%coeffs( i, jp1 )%t_data( 2 ) =          pfac * Eij%coeffs( i, j )%t_data( 1 )
         End Do
      End If

      If( L_max >= 3 ) Then
         Do L = 3, L_max
            Do jp1 = Max( 1, L - 1 - l1 + 1 ), Min( L, l2 )

               j = jp1 - 1
               i = L - 1 - j

               t = 0
               Eij%coeffs( i, jp1 )%t_data( t ) = + qxq_over_a2 * Eij%coeffs( i, j )%t_data( t     ) + &
                                                                  Eij%coeffs( i, j )%t_data( t + 1 )

               t_fac = 1.0_wp
               Do t = 1, L - 2
                  t_fac = t_fac + 1.0_wp
                  Eij%coeffs( i, jp1 )%t_data( t ) =        pfac * Eij%coeffs( i, j )%t_data( t - 1 ) + &
                                                     qxq_over_a2 * Eij%coeffs( i, j )%t_data( t     ) + &
                                                           t_fac * Eij%coeffs( i, j )%t_data( t + 1 )
               End Do

               t = L - 1
               Eij%coeffs( i, jp1 )%t_data( t ) =        pfac * Eij%coeffs( i, j )%t_data( t - 1 ) + &
                                                  qxq_over_a2 * Eij%coeffs( i, j )%t_data( t     )

               t = L
               Eij%coeffs( i, jp1 )%t_data( t ) =        pfac * Eij%coeffs( i, j )%t_data( t - 1 )

            End Do
         End Do
      End If

    End Subroutine hgauss_calc_coeffs_default
    
  End Subroutine hgauss_calc_coeffs

  Pure Function hgauss_calc_product( Eij, l1, l2, x1, x2, x ) Result( gp )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Implicit None

    Real( wp ) :: gp
    
    Class( hgauss_coeffs ), Intent( In ) :: Eij
    Integer               , Intent( In ) :: l1
    Integer               , Intent( In ) :: l2
    Real( wp )            , Intent( In ) :: x1
    Real( wp )            , Intent( In ) :: x2
    Real( wp )            , Intent( In ) :: x

    ! Stupidly large value of L for buffers
    Integer, Parameter :: L_absolute_max = 2000

    Real( wp ), Dimension( 0:L_absolute_max ) :: Hermite
    Real( wp ), Dimension( 0:L_absolute_max ) :: Hermite_prime

    Real( wp ) :: p, xp, pxxp
    Real( wp ) :: sqrt_p
    Real( wp ) :: pfac

    Integer :: L
    Integer :: t

    p  = Eij%p
    xp = ( Eij%a1 * x1 + Eij%a2 * x2 ) / p
    sqrt_p = Sqrt( p )
    pxxp = sqrt_p * ( x - xp )

    L = l1 + l2

    Call Hermite_polynomial( pxxp, Hermite( 0:L ), Hermite_prime( 0:L ) )

    pfac = 1.0_wp
    gp = 0.0_wp
    Do t = 0, L
       gp = gp + Eij%coeffs( l1, l2 )%t_data( t ) * pfac * Hermite( t )
       pfac = pfac * sqrt_p
    End Do

    gp = gp * Exp( - pxxp * pxxp )    
    
  End Function hgauss_calc_product

  Pure Subroutine init_with_NaNs( data )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
    Use, Intrinsic :: ieee_arithmetic, Only : ieee_value, ieee_support_nan, ieee_signaling_nan, &
         ieee_support_halting, ieee_get_halting_mode, ieee_set_halting_mode,                    &
         ieee_get_flag, ieee_set_flag, ieee_invalid

    Implicit None

    Real( wp ), Dimension( : ), Intent( Out ) :: data

    Logical :: is_halting, is_invalid

    If( ieee_support_nan( data ) .And. ieee_support_halting( ieee_invalid ) ) Then

       ! This processor supports ieee maths and control of halting - Use this as carefully as possible to initialise
       ! matrix type objects to a value that can help detect their use when unitilised - namely a ignalling NaN
       ! First get the current halting mode for ieee invalid 
       Call ieee_get_halting_mode( ieee_invalid, is_halting )
       
       ! Now deliberately turn halting off
       Call ieee_set_halting_mode( ieee_invalid, .False. )

       ! Get the current value of the invalid flag to avoid spurious signalling caused by the below
       Call ieee_get_flag( ieee_invalid, is_invalid )

       ! Now can set the data to signalling NaN without worry of stopping the program
       data = ieee_value( data, ieee_signaling_nan )

       ! Reset the invalid flag to what it was before we deliberatley used a signalling NaN to avoid missing ones later
       Call ieee_set_flag( ieee_invalid, is_invalid )

       ! And reset the halting mode
       Call ieee_set_halting_mode( ieee_invalid, is_halting )

    Else

       ! NaNs and/or halting not supported, set to a silly value
       data = Huge( data )

    End If
    
  End Subroutine init_with_NaNs

  Pure Subroutine Hermite_polynomial( x, H, H_prime )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Implicit None

    Real( wp ),                  Intent( In    ) :: x
    Real( wp ), Dimension( 0: ), Intent(   Out ) :: H
    Real( wp ), Dimension( 0: ), Intent(   Out ) :: H_prime
    
    Integer :: n

    H      ( 0 ) = 1.0_wp
    H_prime( 0 ) = 0.0_wp
    Do n = 1, Ubound( H, Dim = 1 )
       H_prime( n ) = 2.0_wp * n * H( n - 1 )
       H      ( n ) = 2.0_wp * x * H( n - 1 ) - H_prime( n - 1 )
    End Do

  End Subroutine Hermite_polynomial
  
End Module hgauss_module
