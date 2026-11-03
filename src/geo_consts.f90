module geo_consts
  use iso_c_binding, only: dp => c_double, i8 => c_int64_t, bool => c_bool
  implicit none
  private
  real(dp), parameter :: PI = 3.141592653589793_dp
  real(dp), parameter :: G = 6.67430e-11_dp ! m3kg-1s-2 tomado de NIST
  public :: G, PI, dp, i8, bool
end module geo_consts
