module geo_utils
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use geo_consts, only: dp, PI
  implicit none
  private

  ! Limites y tolerancias de las latitudes
  real(dp), parameter :: LIM_LAT = 0.5_dp*PI
  real(dp), parameter :: TOL_LAT = 1.0e-12_dp

  public :: deg2rad, rad2deg, validation_lat, lat2colat, colat2lat
contains
  elemental function validation_lat(angle) result(lat_ok)
    real(dp), intent(in) :: angle; real(dp) :: lat_ok
    if (.not. (abs(angle) <= LIM_LAT + TOL_LAT)) then
      lat_ok = ieee_value(angle, ieee_quiet_nan)
    else
      lat_ok = max(-LIM_LAT, min(LIM_LAT, angle))
    end if
  end function validation_lat

  ! Una función elemental funciona igual para escalares y arreglos
  elemental function deg2rad(deg) result(rad)
    real(dp), intent(in) :: deg; real(dp) :: rad
    rad = deg*(PI/180.0_dp)
  end function deg2rad

  elemental function rad2deg(rad) result(deg)
    real(dp), intent(in) :: rad; real(dp) :: deg
    deg = rad*(180.0_dp/PI)
  end function rad2deg

  elemental function lat2colat(lat) result(colat)
    real(dp), intent(in) :: lat; real(dp) :: colat
    colat = 0.5_dp*PI - validation_lat(lat)
  end function lat2colat

  elemental function colat2lat(colat) result(lat)
    real(dp), intent(in) :: colat; real(dp) :: lat
    lat = validation_lat(0.5_dp*PI - colat)
  end function colat2lat
end module geo_utils
