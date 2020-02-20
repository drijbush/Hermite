Module hgauss_module

  ! hermite gaussian module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Private :: array_wrap
     Real( wp ), Dimension( : ), Allocatable, Private :: i_data
  End Type array_wrap

  Type, Private :: array_of_array
     Type( array_wrap ), Dimension( : ), Allocatable, Private :: j_data
  End Type array_of_array

  Type, Public :: hgauss_coeffs
     Private
     Integer                                            , Private :: li_max
     Integer                                            , Private :: lj_max
     Type( array_of_array ), Dimension( : ), Allocatable, Private :: coeffs
     Real( wp )            , Dimension( : ), Allocatable, Private :: H
     Real( wp )            , Dimension( : ), Allocatable, Private :: H_prime
   Contains
     Procedure, Public :: allocate     => hgauss_alloc
     Procedure, Public :: calc_coeffs  => hgauss_calc_coeffs
     Procedure, Public :: calc_product => hgauss_calc_product
  End Type hgauss_coeffs
  
  Private

Contains

  Pure Subroutine hgauss_alloc( Eij, li_max, lj_max, init_with_NaNs )
    
    Class( hgauss_coeffs ), Intent( InOut )           :: Eij
    Integer               , Intent( In    )           :: li_max
    Integer               , Intent( In    )           :: lj_max
    Logical               , Intent( In    ), Optional :: init_with_NaNs

    Call with_dealloc( Eij, li_max, lj_max, init_with_NaNs )

  Contains

    Pure Subroutine with_dealloc( Eij, li_max, lj_max, init_with_NaNs )

      Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
      Use, Intrinsic :: ieee_arithmetic, Only : ieee_value, ieee_support_nan, ieee_signaling_nan, &
           ieee_support_halting, ieee_get_halting_mode, ieee_set_halting_mode,                    &
           ieee_get_flag, ieee_set_flag, ieee_invalid
      Type( hgauss_coeffs ), Intent(   Out )           :: Eij
      Integer              , Intent( In    )           :: li_max
      Integer              , Intent( In    )           :: lj_max
      Logical              , Intent( In    ), Optional :: init_with_NaNs

      Integer :: L
      Integer :: lj, t

      Logical :: is_halting, is_invalid
      Logical :: do_init_with_NaNs

      do_init_with_NaNs = .False.
      If( Present( init_with_NaNs ) ) Then
         If( init_with_NaNs ) Then
            do_init_with_NaNs = ieee_support_nan( 1.0_wp ) .And. ieee_support_halting( ieee_invalid )
         End If
      End If

      Eij%li_max = li_max
      Eij%lj_max = lj_max

      L = li_max + lj_max

      ! Buffers for the values of the hermite polynomials
      Allocate( Eij%H      ( 0:L ) )
      Allocate( Eij%H_prime( 0:L ) )

      Allocate( Eij%coeffs( 0:L ) )
      Do t = 0, L
         Allocate( Eij%coeffs( t )%j_data( 0:lj_max ) )
         Do lj = 0, lj_max
            Allocate( Eij%coeffs( t )%j_data( lj )%i_data( 0:li_max ) )
            If( do_init_with_NaNs ) Then
               ! This processor supports ieee maths and control of halting - Use this as carefully as possible to initialise
               ! matrix type objects to a value that can help detect their use when unitilised - namely a ignalling NaN
               ! First get the current halting mode for ieee invalid 
               Call ieee_get_halting_mode( ieee_invalid, is_halting )
               ! Now deliberately turn halting off
               Call ieee_set_halting_mode( ieee_invalid, .False. )
               ! Get the current value of the invalid flag to avoid spurious signalling caused by the below
               Call ieee_get_flag( ieee_invalid, is_invalid )
               Eij%coeffs( t )%j_data( lj )%i_data = ieee_value( Eij%coeffs( t )%j_data( lj )%i_data, ieee_signaling_nan )
               ! Reset the invalid flag to what it was before we deliberatley used a signalling NaN to avoid missing ones later
               Call ieee_set_flag( ieee_invalid, is_invalid )
               ! And reset the halting mode
               Call ieee_set_halting_mode( ieee_invalid, is_halting )
            End If
         End Do
      End Do

    End Subroutine with_dealloc
    
  End Subroutine hgauss_alloc

  Subroutine hgauss_calc_coeffs( Eij, l1, l2, a1, a2, x1, x2 )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Class( hgauss_coeffs ), Intent( InOut ) :: Eij
    Integer               , Intent( In    ) :: l1
    Integer               , Intent( In    ) :: l2
    Real( wp )            , Intent( In    ) :: a1
    Real( wp )            , Intent( In    ) :: a2
    Real( wp )            , Intent( In    ) :: x1
    Real( wp )            , Intent( In    ) :: x2

    Real( wp ) :: p, q
    Real( wp ) :: xq
    Real( wp ) :: Kab
    Real( wp ) :: pfac, qxq
    Real( wp ) :: qxq_over_a1, qxq_over_a2
    Real( wp ) :: t_fac
   
    Integer :: L_max, L
    Integer :: i, j, t
    Integer :: ip1, jp1
    
    q = a1 * a2 / ( a1 + a2 )

    xq = x1 - x2
    Kab = Exp( - q * xq * xq )

    L_max = l1 + l2

    Eij%coeffs( 0 )%j_data( 0 )%i_data( 0 ) = Kab

    If( L_max > 0 ) Then
   
       p = a1 + a2
       pfac = 0.5_wp / p

       qxq  = q * xq
       qxq_over_a1 = qxq / a1
       qxq_over_a2 = qxq / a2

       ! Fill in the j = 0 values for all values of i
       ! First the unusual low vaues of l1 where some terms are zero since thye just don't exist
       If( l1 >= 1 ) Then
          Eij%coeffs( 0 )%j_data( 0 )%i_data( 1 ) = - qxq_over_a1 * Eij%coeffs( 0 )%j_data( 0 )%i_data( 0 )
          Eij%coeffs( 1 )%j_data( 0 )%i_data( 1 ) =          pfac * Eij%coeffs( 0 )%j_data( 0 )%i_data( 0 )
       End If

       If( l1 >= 2 ) Then
          Eij%coeffs( 0 )%j_data( 0 )%i_data( 2 ) = - qxq_over_a1 * Eij%coeffs( 0 )%j_data( 0 )%i_data( 1 ) + &
                                                                    Eij%coeffs( 1 )%j_data( 0 )%i_data( 1 ) 
          Eij%coeffs( 1 )%j_data( 0 )%i_data( 2 ) =          pfac * Eij%coeffs( 0 )%j_data( 0 )%i_data( 1 ) - &
                                                      qxq_over_a1 * Eij%coeffs( 1 )%j_data( 0 )%i_data( 1 )
          Eij%coeffs( 2 )%j_data( 0 )%i_data( 2 ) =          pfac * Eij%coeffs( 1 )%j_data( 0 )%i_data( 1 )
       End If
       
       ! Now the more general terms
       Do L = 3, l1

          ip1 = L
          i = ip1 - 1

          t_fac = 1.0_wp

          t = 0
          Eij%coeffs( t )%j_data( 0 )%i_data( ip1 ) = - qxq_over_a1 * Eij%coeffs( t     )%j_data( 0 )%i_data( i ) + &
                                                              t_fac * Eij%coeffs( t + 1 )%j_data( 0 )%i_data( i )

          Do t = 1, L - 2
             t_fac = t_fac + 1.0_wp
             Eij%coeffs( t )%j_data( 0 )%i_data( ip1 ) =        pfac * Eij%coeffs( t - 1 )%j_data( 0 )%i_data( i ) - &
                                                         qxq_over_a1 * Eij%coeffs( t     )%j_data( 0 )%i_data( i ) + &
                                                               t_fac * Eij%coeffs( t + 1 )%j_data( 0 )%i_data( i )
          End Do

          t = L - 1
          Eij%coeffs( t )%j_data( 0 )%i_data( ip1 ) =        pfac * Eij%coeffs( t - 1 )%j_data( 0 )%i_data( i ) - &
                                                      qxq_over_a1 * Eij%coeffs( t     )%j_data( 0 )%i_data( i )

          t = L
          Eij%coeffs( t )%j_data( 0 )%i_data( ip1 ) =        pfac * Eij%coeffs( t - 1 )%j_data( 0 )%i_data( i )
          
       End Do

       ! OK All j = 0 values filled in. Now can fill in the rest
       If( l2 >= 1 ) Then
          Eij%coeffs( 0 )%j_data( 1 )%i_data( 0 ) = + qxq_over_a2 * Eij%coeffs( 0 )%j_data( 0 )%i_data( 0 )
          Eij%coeffs( 1 )%j_data( 1 )%i_data( 0 ) =          pfac * Eij%coeffs( 0 )%j_data( 0 )%i_data( 0 )
       End If

       If( L_max >= 2 ) Then
          L = 2
          Do jp1 = Max( 1, L - 1 - l1 + 1 ), Min( L, l2 )
             j = jp1 - 1
             i = L - 1 - j
             Eij%coeffs( 0 )%j_data( jp1 )%i_data( i ) = + qxq_over_a2 * Eij%coeffs( 0 )%j_data( j )%i_data( i ) + &
                                                                         Eij%coeffs( 1 )%j_data( j )%i_data( i )
             Eij%coeffs( 1 )%j_data( jp1 )%i_data( i ) =          pfac * Eij%coeffs( 0 )%j_data( j )%i_data( i ) + &
                                                           qxq_over_a2 * Eij%coeffs( 1 )%j_data( j )%i_data( i )
             Eij%coeffs( 2 )%j_data( jp1 )%i_data( i ) =          pfac * Eij%coeffs( 1 )%j_data( j )%i_data( i )
          End Do
       End If

       If( L_max >= 3 ) Then
          Do L = 3, L_max
             Do jp1 = Max( 1, L - 1 - l1 + 1 ), Min( L, l2 )

                j = jp1 - 1
                i = L - 1 - j

                t = 0
                Eij%coeffs( t )%j_data( jp1 )%i_data( i ) = + qxq_over_a2 * Eij%coeffs( t     )%j_data( j )%i_data( i ) + &
                                                                            Eij%coeffs( t + 1 )%j_data( j )%i_data( i )

                t_fac = 1.0_wp
                Do t = 1, L - 2
                   t_fac = t_fac + 1.0_wp
                   Eij%coeffs( t )%j_data( jp1 )%i_data( i ) =         pfac * Eij%coeffs( t - 1 )%j_data( j )%i_data( i ) + &
                                                                qxq_over_a2 * Eij%coeffs( t     )%j_data( j )%i_data( i ) + &
                                                                      t_fac * Eij%coeffs( t + 1 )%j_data( j )%i_data( i )
                End Do

                t = L - 1
                Eij%coeffs( t )%j_data( jp1 )%i_data( i ) =        pfac * Eij%coeffs( t - 1 )%j_data( j )%i_data( i ) + &
                                                            qxq_over_a2 * Eij%coeffs( t     )%j_data( j )%i_data( i )
                
                t = L
                Eij%coeffs( t )%j_data( jp1 )%i_data( i ) = pfac * Eij%coeffs( t - 1 )%j_data( j )%i_data( i )
                
             End Do
          End Do
       End If

    End If

  End Subroutine hgauss_calc_coeffs

  Function hgauss_calc_product( Eij, l1, l2, a1, a2, x1, x2, x ) Result( gp )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Real( wp ) :: gp
    
    Class( hgauss_coeffs ), Intent( InOut ) :: Eij
    Integer               , Intent( In    ) :: l1
    Integer               , Intent( In    ) :: l2
    Real( wp )            , Intent( In    ) :: a1
    Real( wp )            , Intent( In    ) :: a2
    Real( wp )            , Intent( In    ) :: x1
    Real( wp )            , Intent( In    ) :: x2
    Real( wp )            , Intent( In    ) :: x

    Real( wp ) :: p, xp, pxxp
    Real( wp ) :: sqrt_p
    Real( wp ) :: pfac

    Integer :: L
    Integer :: t

    p = a1 + a2
    xp = ( a1 * x1 + a2 * x2 ) / p
    sqrt_p = Sqrt( p )
    pxxp = sqrt_p * ( x - xp )

    L = l1 + l2

    Call Hermite_polynomial( pxxp, Eij%H( 0:L ), Eij%H_prime( 0:L ) )

    pfac = 1.0_wp
    gp = 0.0_wp
    Do t = 0, L
       gp = gp + Eij%coeffs( t )%j_data( l2 )%i_data( l1 ) * pfac * Eij%H( t )
       pfac = pfac * sqrt_p
    End Do

    gp = gp * Exp( - pxxp * pxxp )
    
    
  End Function hgauss_calc_product

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
