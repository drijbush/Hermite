Program hermite_test

  Use, Intrinsic :: iso_fortran_env, Only : wp => real64
  Use            :: hgauss_module  , Only : hgauss_coeffs
  
  Implicit None

  Integer, Parameter :: l_max = 5

  Type( hgauss_coeffs ) :: HEij
  
  Real( wp ), Dimension( 0:2 * l_max, 0:2 * l_max, 0: 2 * l_max ) :: Eij

  Real( wp ), Dimension( 0:2 * l_max ) :: H, H_prime
  
  Real( wp ) :: a1, a2
  Real( wp ) :: x1, x2
  Real( wp ) :: q, Kab
  Real( wp ) :: p, xp, xq
  Real( wp ) :: g1, g2, gp
  Real( wp ) :: error, max_error
  Real( wp ) :: rel_error, max_rel_error
  Real( wp ) :: x, dx = 0.1_wp
  Real( wp ) :: xl1, xl2

  Integer :: l1, l2, L
  Integer :: t

  Character( Len = * ), Parameter :: format = '( f5.2, 4( 1x, f10.6 ), 2( 1x, g16.8 ) )'

  Write( *, * ) 'l1, l2 ?'
  Read ( *, * ) l1, l2
  Write( *, * ) 'a1, a2 ?'
  Read ( *, * ) a1, a2
  Write( *, * ) 'x1, x2 ?'
  Read ( *, * ) x1, x2

  p = a1 + a2
  q = a1 * a2 / ( a1 + a2 )

  xp = ( a1 * x1 + a2 * x2 ) / p

  xq = x1 - x2
  Kab = Exp( - q * xq * xq )

  Write( *, * ) 'Method 1'
  x = -5.0_wp
  max_error = - Huge( max_error )
  Do While( x <= 5.0_wp )
     xl1 = ( x - x1 ) ** l1
     xl2 = ( x - x2 ) ** l2
     g1 = xl1 * Exp( - a1 * ( x - x1 ) * ( x - x1 ) )
     g2 = xl2 * Exp( - a2 * ( x - x2 ) * ( x - x2 ) )
     gp = xl1 * xl2 * Kab * Exp( - p * ( x - xp ) * ( x - xp ) )
     error = Abs( g1 * g2 - gp )
     max_error = Max( error, max_error )
     rel_error = error / Abs( g1 * g2 )
     max_rel_error = Max( rel_error, max_rel_error )
     Write( *, format ) x, g1, g2, g1 * g2, gp, error, rel_error
     x = x + dx
  End Do
  Write( *, * ) max_error, max_rel_error

  Write( *, * ) 'Method 2'
  L = l1 + l2
  Call calc_hermite_expansion_coeffs( l1, l2, a1, a2, x1, x2, Eij )
  x = -5.0_wp
  max_error = - Huge( max_error )
  Do While( x <= 2.8_wp )
     xl1 = ( x - x1 ) ** l1
     xl2 = ( x - x2 ) ** l2
     g1 = xl1 * Exp( - a1 * ( x - x1 ) * ( x - x1 ) )
     g2 = xl2 * Exp( - a2 * ( x - x2 ) * ( x - x2 ) )
     gp = 0.0_wp
     Call Hermite_polynomial( Sqrt( p ) * ( x - xp ), H( 0:L ), H_prime( 0:L ) )
     Do t = 0, l1 + l2
        gp = gp + &
             Eij( l1, l2, t ) * ( p ** ( 0.5_wp * t ) ) * H( t )
     End Do
     gp = gp * Exp( - p * ( x - xp ) * ( x - xp ) )
     error = Abs( g1 * g2 - gp )
     max_error = Max( error, max_error )
     rel_error = error / Abs( g1 * g2 )
     max_rel_error = Max( rel_error, max_rel_error )
     Write( *, format ) x, g1, g2, g1 * g2, gp, error, rel_error
     x = x + dx
  End Do
  Write( *, * ) max_error, max_rel_error

  Write( *, '( 100( g12.4, 1x ) )' ) Eij( l1, l2, : )

  Write( *, * ) 'Method 3'
  Call HEij%allocate( l1, l2, .True. )
  Call HEij%calc_coeffs( l1, l2, a1, a2, x1, x2 )
  L = l1 + l2
  x = -5.0_wp
  max_error = - Huge( max_error )
  Do While( x <= 5.0_wp )
     xl1 = ( x - x1 ) ** l1
     xl2 = ( x - x2 ) ** l2
     g1 = xl1 * Exp( - a1 * ( x - x1 ) * ( x - x1 ) )
     g2 = xl2 * Exp( - a2 * ( x - x2 ) * ( x - x2 ) )
     gp = HEij%calc_product( l1, l2, a1, a2, x1, x2, x )
     error = Abs( g1 * g2 - gp )
     max_error = Max( error, max_error )
     rel_error = error / Abs( g1 * g2 )
     max_rel_error = Max( rel_error, max_rel_error )
     Write( *, format ) x, g1, g2, g1 * g2, gp, error, rel_error
     x = x + dx
  End Do
  Write( *, * ) max_error, max_rel_error

Contains

  Pure Subroutine calc_hermite_expansion_coeffs( l1, l2, a1, a2, x1, x2, Eij )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
    Use, Intrinsic :: ieee_arithmetic, Only : ieee_value, ieee_quiet_nan
    
    Implicit None

    Integer                            , Intent( In    ) :: l1
    Integer                            , Intent( In    ) :: l2
    Real( wp )                         , Intent( In    ) :: a1
    Real( wp )                         , Intent( In    ) :: a2
    Real( wp )                         , Intent( In    ) :: x1
    Real( wp )                         , Intent( In    ) :: x2
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent(   Out ) :: Eij

    Real( wp ) :: p, q
    Real( wp ) :: xq
    Real( wp ) :: Kab

    Integer :: L_max, L
    Integer :: i, j, t
    Integer :: ip1, jp1
    
    p = a1 + a2
    q = a1 * a2 / ( a1 + a2 )

    xq = x1 - x2
    Kab = Exp( - q * xq * xq )

!!$    Eij( 0:l1, 0:l2, 0:l1 + l2 ) = Huge( Eij )
    Eij = ieee_value( Eij, ieee_quiet_nan )

    ! Start off manually
    ! L = l1 + l2 = 0
    ! t = 0
    Eij( 0, 0, 0 ) = Kab

    ! L = 1
    ! t = 0
    Eij( 1, 0, 0 ) = - q * xq * Eij( 0, 0, 0 ) / a1
    Eij( 0, 1, 0 ) = + q * xq * Eij( 0, 0, 0 ) / a2

    ! t = 1
    Eij( 1, 0, 1 ) = Eij( 0, 0, 0 ) / ( 2.0_wp * p )
    Eij( 0, 1, 1 ) = Eij( 0, 0, 0 ) / ( 2.0_wp * p )

    ! L = 2
    ! t = 0
    Eij( 2, 0, 0 ) = - q * xq * Eij( 1, 0, 0 ) / a1 + 1.0_wp * Eij( 1, 0, 1 )
    Eij( 1, 1, 0 ) = + q * xq * Eij( 1, 0, 0 ) / a2 + 1.0_wp * Eij( 1, 0, 1 )
    Eij( 0, 2, 0 ) = + q * xq * Eij( 0, 1, 0 ) / a2 + 1.0_wp * Eij( 0, 1, 1 )
    
    ! t = 1
    Eij( 2, 0, 1 ) = Eij( 1, 0, 0 ) / ( 2.0_wp * p ) - q * xq * Eij( 1, 0, 1 ) / a1
    Eij( 1, 1, 1 ) = Eij( 1, 0, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 1, 0, 1 ) / a2
    Eij( 0, 2, 1 ) = Eij( 0, 1, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 0, 1, 1 ) / a2

    ! t = 2
    Eij( 2, 0, 2 ) = Eij( 1, 0, 1 ) / ( 2.0_wp * p )
    Eij( 1, 1, 2 ) = Eij( 1, 0, 1 ) / ( 2.0_wp * p )
    Eij( 0, 2, 2 ) = Eij( 0, 1, 1 ) / ( 2.0_wp * p )

    ! L = 3
    ! t = 0
    Eij( 3, 0, 0 ) = - q * xq * Eij( 2, 0, 0 ) / a1 + 1.0_wp * Eij( 2, 0, 1 )
    Eij( 2, 1, 0 ) = + q * xq * Eij( 2, 0, 0 ) / a2 + 1.0_wp * Eij( 2, 0, 1 )
    Eij( 1, 2, 0 ) = + q * xq * Eij( 1, 1, 0 ) / a2 + 1.0_wp * Eij( 1, 1, 1 )
    Eij( 0, 3, 0 ) = + q * xq * Eij( 0, 2, 0 ) / a2 + 1.0_wp * Eij( 0, 2, 1 )

    ! t = 1
    Eij( 3, 0, 1 ) = Eij( 2, 0, 0 ) / ( 2.0_wp * p ) - q * xq * Eij( 2, 0, 1 ) / a1 + 2.0_wp * Eij( 2, 0, 2 )
    Eij( 2, 1, 1 ) = Eij( 2, 0, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 2, 0, 1 ) / a2 + 2.0_wp * Eij( 2, 0, 2 )
    Eij( 1, 2, 1 ) = Eij( 1, 1, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 1, 1, 1 ) / a2 + 2.0_wp * Eij( 1, 1, 2 )
    Eij( 0, 3, 1 ) = Eij( 0, 2, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 0, 2, 1 ) / a2 + 2.0_wp * Eij( 0, 2, 2 )

    ! t = 2
    Eij( 3, 0, 2 ) = Eij( 2, 0, 1 ) / ( 2.0_wp * p ) - q * xq * Eij( 2, 0, 2 ) / a1 
    Eij( 2, 1, 2 ) = Eij( 2, 0, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 2, 0, 2 ) / a2 
    Eij( 1, 2, 2 ) = Eij( 1, 1, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 1, 1, 2 ) / a2 
    Eij( 0, 3, 2 ) = Eij( 0, 2, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 0, 2, 2 ) / a2

    ! t = 3
    Eij( 3, 0, 3 ) = Eij( 2, 0, 2 ) / ( 2.0_wp * p )
    Eij( 2, 1, 3 ) = Eij( 2, 0, 2 ) / ( 2.0_wp * p )
    Eij( 1, 2, 3 ) = Eij( 1, 1, 2 ) / ( 2.0_wp * p )
    Eij( 0, 3, 3 ) = Eij( 0, 2, 2 ) / ( 2.0_wp * p )

    ! L = 4
    ! t = 0
    Eij( 4, 0, 0 ) = - q * xq * Eij( 3, 0, 0 ) / a1 + 1.0_wp * Eij( 3, 0, 1 )
    Eij( 3, 1, 0 ) = + q * xq * Eij( 3, 0, 0 ) / a2 + 1.0_wp * Eij( 3, 0, 1 )
    Eij( 2, 2, 0 ) = + q * xq * Eij( 2, 1, 0 ) / a2 + 1.0_wp * Eij( 2, 1, 1 )
    Eij( 1, 3, 0 ) = + q * xq * Eij( 1, 2, 0 ) / a2 + 1.0_wp * Eij( 1, 2, 1 )
    Eij( 0, 4, 0 ) = + q * xq * Eij( 0, 3, 0 ) / a2 + 1.0_wp * Eij( 0, 3, 1 )

    ! t = 1
    Eij( 4, 0, 1 ) = Eij( 3, 0, 0 ) / ( 2.0_wp * p ) - q * xq * Eij( 3, 0, 1 ) / a1 + 2.0_wp * Eij( 3, 0, 2 )
    Eij( 3, 1, 1 ) = Eij( 3, 0, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 3, 0, 1 ) / a2 + 2.0_wp * Eij( 3, 0, 2 )
    Eij( 2, 2, 1 ) = Eij( 2, 1, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 2, 1, 1 ) / a2 + 2.0_wp * Eij( 2, 1, 2 )
    Eij( 1, 3, 1 ) = Eij( 1, 2, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 1, 2, 1 ) / a2 + 2.0_wp * Eij( 1, 2, 2 )
    Eij( 0, 4, 1 ) = Eij( 0, 3, 0 ) / ( 2.0_wp * p ) + q * xq * Eij( 0, 3, 1 ) / a2 + 2.0_wp * Eij( 0, 3, 2 )
    
    ! t = 2
    Eij( 4, 0, 2 ) = Eij( 3, 0, 1 ) / ( 2.0_wp * p ) - q * xq * Eij( 3, 0, 2 ) / a1 + 3.0_wp * Eij( 3, 0, 3 )
    Eij( 3, 1, 2 ) = Eij( 3, 0, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 3, 0, 2 ) / a2 + 3.0_wp * Eij( 3, 0, 3 )
    Eij( 2, 2, 2 ) = Eij( 2, 1, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 2, 1, 2 ) / a2 + 3.0_wp * Eij( 2, 1, 3 )
    Eij( 1, 3, 2 ) = Eij( 1, 2, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 1, 2, 2 ) / a2 + 3.0_wp * Eij( 1, 2, 3 )
    Eij( 0, 4, 2 ) = Eij( 0, 3, 1 ) / ( 2.0_wp * p ) + q * xq * Eij( 0, 3, 2 ) / a2 + 3.0_wp * Eij( 0, 3, 3 )

    ! t = 3
    Eij( 4, 0, 3 ) = Eij( 3, 0, 2 ) / ( 2.0_wp * p ) - q * xq * Eij( 3, 0, 3 ) / a1 
    Eij( 3, 1, 3 ) = Eij( 3, 0, 2 ) / ( 2.0_wp * p ) + q * xq * Eij( 3, 0, 3 ) / a2 
    Eij( 2, 2, 3 ) = Eij( 2, 1, 2 ) / ( 2.0_wp * p ) + q * xq * Eij( 2, 1, 3 ) / a2 
    Eij( 1, 3, 3 ) = Eij( 1, 2, 2 ) / ( 2.0_wp * p ) + q * xq * Eij( 1, 2, 3 ) / a2
    Eij( 0, 4, 3 ) = Eij( 0, 3, 2 ) / ( 2.0_wp * p ) + q * xq * Eij( 0, 3, 3 ) / a2
    
    ! t = 4
    Eij( 4, 0, 4 ) = Eij( 3, 0, 3 ) / ( 2.0_wp * p )
    Eij( 3, 1, 4 ) = Eij( 3, 0, 3 ) / ( 2.0_wp * p )
    Eij( 2, 2, 4 ) = Eij( 2, 1, 3 ) / ( 2.0_wp * p )
    Eij( 1, 3, 4 ) = Eij( 1, 2, 3 ) / ( 2.0_wp * p )
    Eij( 0, 4, 4 ) = Eij( 0, 3, 3 ) / ( 2.0_wp * p )

!!$    l_max = Ubound( Eij, Dim = 1 )
    l_max = l1 + l2
    Do L = 5, l_max

       ! t = 0
       t = 0
       ip1 = L
       i = ip1 - 1
       Eij( ip1, 0, t ) = - q * xq * Eij( i, 0, t ) / a1 + ( t + 1 ) * Eij(  i, 0, t + 1 )
       Do jp1 = 1, L
          j = jp1 - 1
          i = L - 1 - j
          Eij( i, jp1, t ) = + q * xq * Eij( i, j, t ) / a2 + ( t + 1 ) * Eij( i, j, t + 1 )
       End Do

       ! t = 1 ... L - 2
       Do t = 1, L - 2
          ip1 = L
          i = ip1 - 1
          Eij( ip1, 0, t ) = Eij( i, 0, t - 1 ) / ( 2.0_wp * p ) - q * xq * Eij( i, 0, t ) / a1 + ( t + 1 ) * Eij( i, 0, t + 1 )
          Do jp1 = 1, L
             j = jp1 - 1
             i = L - 1 - j
             Eij( i, jp1, t ) = Eij( i, j, t - 1 ) / ( 2.0_wp * p ) + q * xq * Eij( i, j, t ) / a2 + ( t + 1 ) * Eij( i, j, t + 1 )
          End Do
       End Do

       ! t = L - 1
       t = L - 1
       ip1 = L
       i = ip1 - 1
       Eij( ip1, 0, t ) = Eij( i, 0, t - 1 ) / ( 2.0_wp * p ) - q * xq * Eij( i, 0, t ) / a1 
       Do jp1 = 1, L
          j = jp1 - 1
          i = L - 1 - j
          Eij( i, jp1, t ) = Eij( i, j, t - 1 ) / ( 2.0_wp * p ) + q * xq * Eij( i, j, t ) / a2
       End Do

       ! t = L
       t = L
       ip1 = L
       i = ip1 - 1
       Eij( ip1, 0, t ) = Eij( i, 0, t - 1 ) / ( 2.0_wp * p )
       Do jp1 = 1, L
          j = jp1 - 1
          i = L - 1 - j
          Eij( i, jp1, t ) = Eij( i, j, t - 1 ) / ( 2.0_wp * p )
       End Do
       
    End Do

  End Subroutine calc_hermite_expansion_coeffs

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

!!$    Do n = 1, Ubound( H, Dim = 1 )
!!$       H_prime( n ) = n * H( n - 1 )
!!$       H      ( n ) = x * H( n - 1 ) - H_prime( n - 1 ) * 0.5_wp
!!$    End Do

!!$    ! Rescale to proper values
!!$    Do n = 1, Ubound( H, Dim = 1 )
!!$       H_prime( n ) = H_prime( n ) * ( 2 ** n )
!!$       H( n )       = H      ( n ) * ( 2 ** n )
!!$    End Do
    
  End Subroutine Hermite_polynomial
  
End Program hermite_test
