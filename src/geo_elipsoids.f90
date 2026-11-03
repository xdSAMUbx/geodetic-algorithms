module geo_elipsoids
  use quadpack, only: dqags
  use geo_utils, only: validation_lat
  use geo_consts, only: dp, PI
  implicit none
  private
  type elipsoid
    real(dp) :: a      ! Semieje mayor [m]
    real(dp) :: inv_f  ! Inverso de aplanamiento
  contains
    ! Obtención de constantes geométricas según el elipsoide
    procedure :: f => achatamiento
    procedure :: b => semieje_menor
    procedure :: ecc => primera_excentricidad      ! e^2
    procedure :: E => excentricidad_lineal
    procedure :: c => radio_polar_curvatura
    procedure :: ecc_2 => segunda_excentricidad    ! e'^2

    ! Funciones auxiliares dependientes de la latitud geodésica
    procedure :: W => w_geodesica                  ! W = sqrt(1 - e^2 sin^2(phi))
    procedure :: V => v_geodesica                  ! V = sqrt(1 + e'^2 cos^2(phi))

    ! Conversión de Latitudes
    procedure :: geod2red
    procedure :: red2geod
    procedure :: geod2geoc
    procedure :: geoc2geod
    procedure :: red2geoc
    procedure :: geoc2red

    ! Calculo de Radios de Curvatura
    procedure :: M => radio_meridiano          ! Radio de curvatura en el meridiano
    procedure :: N => radio_primer_vertical    ! Radio de curvatura vertical primario

    ! Cuadrante del meridiano y radios medios de la Tierra
    procedure :: Q => cuadrante_meridiano      ! c * Int dphi/V^3
    procedure :: R1                            ! Radio medio aritmético
    procedure :: R2                            ! Radio de la esfera de igual área
    procedure :: R3                            ! Radio de la esfera de igual volumen
  end type elipsoid

  type(elipsoid), parameter :: GRS80 = elipsoid(6378137.0_dp, 298.257222101_dp)
  type(elipsoid), parameter :: WGS84 = elipsoid(6378137.0_dp, 298.257223563_dp)

  ! Elipsoide "activo" para los integrandos: dqags solo llama a f(x), así que
  ! cada rutina guarda aquí una copia del elipsoide antes de integrar.
  ! (No es thread-safe: dos hilos integrando a la vez se pisarían.)
  type(elipsoid), save :: elip_int

  public :: GRS80, WGS84, elipsoid
contains
  pure function achatamiento(self) result(f)
    class(elipsoid), intent(in) :: self
    real(dp) :: f
    f = 1.0_dp/self%inv_f
  end function achatamiento

  pure function semieje_menor(self) result(b)
    class(elipsoid), intent(in) :: self
    real(dp) :: b
    b = self%a*(1.0_dp - self%f())
  end function semieje_menor

  ! Primera excentricidad AL CUADRADO: e^2
  pure function primera_excentricidad(self) result(ecc)
    class(elipsoid), intent(in) :: self
    real(dp) :: ecc; real(dp) :: f
    f = self%f()
    ecc = f*(2.0_dp - f)
  end function primera_excentricidad

  pure function excentricidad_lineal(self) result(E)
    class(elipsoid), intent(in) :: self
    real(dp) :: E; real(dp) :: ecc
    ecc = self%ecc()
    E = self%a*sqrt(ecc)
  end function excentricidad_lineal

  ! Radio polar de curvatura: c = a^2/b
  pure function radio_polar_curvatura(self) result(c)
    class(elipsoid), intent(in) :: self
    real(dp) :: c; real(dp) :: b
    b = self%b()
    c = (self%a**2)/b
  end function radio_polar_curvatura

  ! Segunda excentricidad AL CUADRADO: e'^2 = e^2/(1 - e^2)
  pure function segunda_excentricidad(self) result(ecc_2)
    class(elipsoid), intent(in) :: self
    real(dp) :: ecc_2; real(dp) :: ecc
    ecc = self%ecc()
    ecc_2 = ecc/(1.0_dp - ecc)
  end function segunda_excentricidad

  ! W = sqrt(1 - e^2 sin^2(phi))   (latitud geodésica en radianes)
  elemental function w_geodesica(self, lat) result(w)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: lat
    real(dp) :: w
    w = sqrt(1.0_dp - self%ecc()*sin(validation_lat(lat))**2)
  end function w_geodesica

  ! V = sqrt(1 + e'^2 cos^2(phi))  (latitud geodésica en radianes)
  elemental function v_geodesica(self, lat) result(v)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: lat
    real(dp) :: v
    v = sqrt(1.0_dp + self%ecc_2()*cos(validation_lat(lat))**2)
  end function v_geodesica

  ! Conversión de latitudes
  ! Geodesica -> Reducida
  elemental function geod2red(self, geodesica) result(reducida)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: geodesica
    real(dp) :: reducida, phi

    phi = validation_lat(geodesica)
    reducida = atan2(sqrt(1.0_dp - self%ecc())*sin(phi), cos(phi))
  end function geod2red

  ! Reducida -> Geodesica
  elemental function red2geod(self, reducida) result(geodesica)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: reducida; real(dp) :: geodesica, red

    red = validation_lat(reducida)
    geodesica = atan2(sin(red), cos(red)*sqrt(1.0_dp - self%ecc()))
  end function red2geod

  ! Geodesica -> Geocentrica
  elemental function geod2geoc(self, geodesica) result(geocentrica)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: geodesica; real(dp) :: geocentrica, geod

    geod = validation_lat(geodesica)
    geocentrica = atan2(sin(geod)*(1.0_dp - self%ecc()), cos(geod))
  end function geod2geoc

  ! Geocentrica -> Geodesica
  elemental function geoc2geod(self, geocentrica) result(geodesica)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: geocentrica; real(dp) :: geodesica, geoc

    geoc = validation_lat(geocentrica)
    geodesica = atan2(sin(geoc), (1.0_dp - self%ecc())*cos(geoc))
  end function geoc2geod

  ! Reducida -> Geocentrica
  elemental function red2geoc(self, reducida) result(geocentrica)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: reducida; real(dp) :: geocentrica, red

    red = validation_lat(reducida)
    geocentrica = atan2(sqrt(1.0_dp - self%ecc())*sin(red), cos(red))
  end function red2geoc

  ! Geocentrica -> Reducida
  elemental function geoc2red(self, geocentrica) result(reducida)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: geocentrica; real(dp) :: reducida, geoc

    geoc = validation_lat(geocentrica)
    reducida = atan2(sin(geoc), sqrt(1.0_dp - self%ecc())*cos(geoc))
  end function geoc2red

  ! Radio de curvatura en el meridiano: M = a(1 - e^2)/W^3
  elemental function radio_meridiano(self, geodesica) result(radio)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: geodesica
    real(dp) :: radio

    radio = self%a*(1.0_dp - self%ecc())/self%W(geodesica)**3
  end function radio_meridiano

  ! Radio de curvatura en el primer vertical: N = a/W
  elemental function radio_primer_vertical(self, geodesica) result(radio)
    class(elipsoid), intent(in) :: self
    real(dp), intent(in) :: geodesica
    real(dp) :: radio

    radio = self%a/self%W(geodesica)
  end function radio_primer_vertical

  ! ---------------------------------------------------------------------------
  ! Integrandos (dqags solo sabe llamar a f(x); el elipsoide viaja en elip_int)
  ! ---------------------------------------------------------------------------

  ! 1/V^3 = 1/(1 + e'^2 cos^2(phi))^(3/2)
  function integrando_Q(phi) result(g)
    real(dp), intent(in) :: phi
    real(dp) :: g
    g = 1.0_dp/elip_int%V(phi)**3
  end function integrando_Q

  ! cos(phi)/V^4 = cos(phi)/(1 + e'^2 cos^2(phi))^2
  function integrando_R2(phi) result(g)
    real(dp), intent(in) :: phi
    real(dp) :: g
    g = cos(phi)/elip_int%V(phi)**4
  end function integrando_R2

  ! Arco de meridiano:  Q = c * Integral_{lat1}^{lat2} dphi / V^3   [m]
  ! Forma estándar (sin lat1/lat2): cuadrante, de 0 a pi/2.
  ! Con lat1 > lat2 el resultado sale negativo (arco con signo).
  subroutine cuadrante_meridiano(self, s, abserr, ier, lat1, lat2, epsrel)
    class(elipsoid), intent(in) :: self
    real(dp), intent(out) :: s, abserr
    integer, intent(out) :: ier
    real(dp), intent(in), optional :: lat1, lat2, epsrel

    integer, parameter :: limit = 100, lenw = 4*limit
    integer :: neval, last, iwork(limit)
    real(dp) :: work(lenw), tol, integral, lo, hi

    lo = 0.0_dp; if (present(lat1)) lo = lat1
    hi = 0.5_dp*PI; if (present(lat2)) hi = lat2
    tol = 1.0e-13_dp; if (present(epsrel)) tol = epsrel

    elip_int = self

    call dqags(integrando_Q, lo, hi, 0.0_dp, tol, integral, abserr, &
               neval, ier, limit, lenw, last, iwork, work)

    s = self%c()*integral
    abserr = self%c()*abserr   ! dqags integró 1/V^3 sin el factor c
  end subroutine cuadrante_meridiano

  ! Radio medio (media aritmética de los semiejes): R1 = (2a + b)/3 = a(1 - f/3)
  elemental function R1(self) result(val)
    class(elipsoid), intent(in) :: self
    real(dp) :: val
    val = self%a*(1.0_dp - self%f()/3.0_dp)
  end function R1

  ! R2 = c * sqrt( Integral_{lat1}^{lat2} cos(phi)/V^4 dphi )
  ! Forma estándar (sin lat1/lat2): de 0 a pi/2 -> radio de la esfera de igual área.
  ! Con otros límites ya no es el radio de igual área del elipsoide completo:
  ! es R2 = sqrt(A_zona/(2*pi)), con A_zona el área de la zona entre ambas latitudes.
  subroutine R2(self, r, abserr, ier, lat1, lat2, epsrel)
    class(elipsoid), intent(in) :: self
    real(dp), intent(out) :: r, abserr
    integer, intent(out) :: ier
    real(dp), intent(in), optional :: lat1, lat2, epsrel

    integer, parameter :: limit = 100, lenw = 4*limit
    integer :: neval, last, iwork(limit)
    real(dp) :: work(lenw), tol, integral, err_int, lo, hi, tmp

    lo = 0.0_dp; if (present(lat1)) lo = lat1
    hi = 0.5_dp*PI; if (present(lat2)) hi = lat2
    tol = 1.0e-13_dp; if (present(epsrel)) tol = epsrel

    ! La raíz no distingue el orden de los límites
    if (lo > hi) then
      tmp = lo; lo = hi; hi = tmp
    end if

    ! Intervalo vacío: evita dividir por cero al propagar el error
    if (lo == hi) then
      r = 0.0_dp; abserr = 0.0_dp; ier = 0
      return
    end if

    elip_int = self

    call dqags(integrando_R2, lo, hi, 0.0_dp, tol, integral, err_int, &
               neval, ier, limit, lenw, last, iwork, work)

    r = self%c()*sqrt(integral)
    abserr = self%c()*err_int/(2.0_dp*sqrt(integral))   ! error propagado por la raíz
  end subroutine R2

  ! Radio de la esfera de igual volumen: R3 = (a^2 b)^(1/3)
  elemental function R3(self) result(val)
    class(elipsoid), intent(in) :: self
    real(dp) :: val
    val = (self%a**2*self%b())**(1.0_dp/3.0_dp)
  end function R3
end module geo_elipsoids
