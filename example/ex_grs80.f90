program main
  use geo_elipsoids, only: GRS80
  use geo_consts, only: dp, PI
  implicit none

  real(dp), parameter :: RAD2DEG = 180.0_dp/PI
  real(dp), parameter :: PHI = PI/4.0_dp      ! latitud de ejemplo: 45 grados
  real(dp) :: lat(7), q, q_err, r2, r2_err
  integer :: i, ier_q, ier_r2

  ! Latitudes de la tabla: 0, 15, 30, 45, 60, 75 y 90 grados
  lat(1:6) = [(real(15*(i - 1), dp), i=1, 6)]*PI/180.0_dp
  lat(7) = 0.5_dp*PI

  ! Las dos integrales numéricas (QUADPACK, dqags)
  call GRS80%Q(q, q_err, ier_q)
  call GRS80%R2(r2, r2_err, ier_r2)

  print '(a)', '=========================================================='
  print '(a)', '  ELIPSOIDE GRS80'
  print '(a)', '=========================================================='

  print '(/,a)', '--- Parámetros definidores ---'
  print '(a,f18.6,a)', '  a  (semieje mayor)         = ', GRS80%a, ' m'
  print '(a,f18.9)', '  1/f (inverso aplanamiento) = ', GRS80%inv_f

  print '(/,a)', '--- Constantes geométricas ---'
  print '(a,es22.14)', '  f   (aplanamiento)         = ', GRS80%f()
  print '(a,f18.6,a)', '  b   (semieje menor)        = ', GRS80%b(), ' m'
  print '(a,es22.14)', '  e^2  (1ra excentricidad^2) = ', GRS80%ecc()
  print '(a,es22.14)', '  e''^2 (2da excentricidad^2) = ', GRS80%ecc_2()
  print '(a,f18.6,a)', '  E   (excentricidad lineal) = ', GRS80%E(), ' m'
  print '(a,f18.6,a)', '  c   (radio polar curvat.)  = ', GRS80%c(), ' m'

  print '(/,a)', '--- Cuadrante del meridiano y radios medios ---'
  print '(a,f18.6,a,es9.2,a,i0)', '  Q  (cuadrante)             = ', q, ' m   err=', q_err, '  ier=', ier_q
  print '(a,f18.6,a)', '  R1 (radio medio)           = ', GRS80%R1(), ' m'
  print '(a,f18.6,a,es9.2,a,i0)', '  R2 (igual área)            = ', r2, ' m   err=', r2_err, '  ier=', ier_r2
  print '(a,f18.6,a)', '  R3 (igual volumen)         = ', GRS80%R3(), ' m'

  print '(/,a)', '--- Funciones auxiliares en latitud 45 grados ---'
  print '(a,f18.12)', '  W  = sqrt(1 - e^2 sin^2)   = ', GRS80%W(PHI)
  print '(a,f18.12)', '  V  = sqrt(1 + e''^2 cos^2)  = ', GRS80%V(PHI)
  print '(a,f18.6,a)', '  M  (radio meridiano)       = ', GRS80%M(PHI), ' m'
  print '(a,f18.6,a)', '  N  (radio primer vertical) = ', GRS80%N(PHI), ' m'

  print '(/,a)', '--- Latitudes equivalentes a geodésica = 45 grados ---'
 print '(a,f18.10,a)', '  Geodésica                  = ', PHI*RAD2DEG, ' grados'
  print '(a,f18.10,a)', '  Reducida                   = ', GRS80%geod2red(PHI)*RAD2DEG, ' grados'
  print '(a,f18.10,a)', '  Geocéntrica                = ', GRS80%geod2geoc(PHI)*RAD2DEG, ' grados'

  print '(/,a)', '--- Tabla por latitud geodésica ---'
  print '(a)', '   lat [deg]        W            V               M [m]           N [m]'
  print '(a)', '  ----------------------------------------------------------------------'
  do i = 1, size(lat)
    print '(2x,f8.2,2x,f11.8,2x,f11.8,2x,f15.4,2x,f15.4)', &
          lat(i)*RAD2DEG, GRS80%W(lat(i)), GRS80%V(lat(i)), GRS80%M(lat(i)), GRS80%N(lat(i))
  end do
  print '(a)', '=========================================================='
end program main
