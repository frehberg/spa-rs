// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0>.
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Solar Position Algorithm module for Rust
//!
//! Collection of algorithms calculating sunrise/sunset and azimuth/zenith-angle.
//! The SPA library supports `std` and `no_std` build targets. The platform specific floating operations
//! are to be defined as trait impl `FloatOps`, the example below illustrates the usage for the `std` target build.
//!
//! ```rust
//!
//! use chrono::{TimeZone, Utc};
//! use spa::{solar_position, sunrise_and_set, SolarPos, FloatOps, SunriseAndSet};
//!
//!
//! // FloatOps for the std environment
//! pub struct StdFloatOps;
//!
//! // FloatOps for the std environment, mapping directly onto f64 operations
//! impl FloatOps for StdFloatOps {
//!     fn sin(x: f64) -> f64 { x.sin() }
//!     fn cos(x: f64) -> f64 { x.cos() }
//!     fn tan(x: f64) -> f64 { x.tan() }
//!     fn asin(x: f64) -> f64 { x.asin() }
//!     fn acos(x: f64) -> f64 { x.acos() }
//!     fn atan(x: f64) -> f64 { x.atan() }
//!     fn atan2(y: f64, x: f64) -> f64 { y.atan2(x) }
//!     fn trunc(x: f64) -> f64 { x.trunc() }
//! }
//!
//! // main
//! fn main() {
//!     let datetime = Utc.with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
//!         .single().unwrap();
//!
//!     // geo-pos near Frankfurt/Germany
//!     let lat = 50.0;
//!     let lon = 10.0;
//!
//!     let _solpos: SolarPos = solar_position::<StdFloatOps>(dt, lat, lon).unwrap();
//!     // ...
//!     let _sunrise_set: SunriseAndSet =  sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap();
//!     // ...
//! }
//! ```

#![cfg_attr(not(test), no_std)]

use chrono::prelude::Utc;
use chrono::DateTime;
use chrono::TimeZone;
use chrono::Timelike;
use core::f64::consts::PI;

const PI2: f64 = PI * 2.0;
const RAD: f64 = 0.017453292519943295769236907684886;
const EARTH_MEAN_RADIUS: f64 = 6371.01;
// In km
const ASTRONOMICAL_UNIT: f64 = 149597890.0;
// In km
const JD2000: f64 = 2451545.0;

/// platform specific floating operations
///
/// For example implement as
/// ```rust
/// /// FloatOps type
/// pub struct StdFloatOps;
///
/// /// FloatOps for the std environment, mapping directly onto f64 operations
/// impl FloatOps for StdFloatOps {
///     fn sin(x: f64) -> f64 { x.sin() }
///     fn cos(x: f64) -> f64 { x.cos() }
///     fn tan(x: f64) -> f64 { x.tan() }
///     fn asin(x: f64) -> f64 { x.asin() }
///     fn acos(x: f64) -> f64 { x.acos() }
///     fn atan(x: f64) -> f64 { x.atan() }
///     fn atan2(y: f64, x: f64) -> f64 { y.atan2(x) }
///     fn trunc(x: f64) -> f64 { x.trunc() }
/// }
/// ```
pub trait FloatOps {
    fn sin(x: f64) -> f64;
    fn cos(x: f64) -> f64;
    fn tan(x: f64) -> f64;

    fn asin(x: f64) -> f64;
    fn acos(x: f64) -> f64;
    fn atan(x: f64) -> f64;
    fn atan2(y: f64, x: f64) -> f64;

    fn trunc(x: f64) -> f64;
}


/// The sun-rise and sun-set as UTC, otherwise permanent polar-night or polar-day
#[derive(Debug, Clone)]
pub enum SunriseAndSet {
    /// The polar night occurs in the northernmost and southernmost regions of the Earth when
    /// the night lasts for more than 24 hours. This occurs only inside the polar circles.
    PolarNight,
    /// The polar day occurs when the Sun stays above the horizon for more than 24 hours.
    /// This occurs only inside the polar circles.
    PolarDay,
    /// The sunrise and sunset as UTC, the Coordinated Universal Time (offset +00:00), incl.
    /// leap-seconds.
    ///
    /// These UTC time-points can be transformed to local times based on longitude or the
    /// corresponding timezone. As timezones are corresponding to borders of countries, a map
    /// would be required.
    Daylight(DateTime<Utc>, DateTime<Utc>),
}

/// The solar position
#[derive(Debug, Clone)]
pub struct SolarPos {
    /// horizontal angle measured clockwise from a north base line or meridian
    pub azimuth: f64,
    /// the angle between the zenith and the centre of the sun's disc
    pub zenith_angle: f64,
}

/// The error conditions
#[derive(Debug, Clone)]
pub enum SpaError {
    BadParam,
}

/// Displaying SpaError, enabling error message on console
impl core::fmt::Display for SpaError {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match *self {
            SpaError::BadParam => write!(f, "Latitude or longitude are not within valid ranges."),
        }
    }
}

/// Converting DateTime<Utc> to Julian-Days (f64)
///
/// Julian-Days is the number of days (incl. fraction) since January 1, 4713 BC, noon, without the
/// leap seconds.
///
/// # Arguments
///
/// * `utc` - UTC time point
///
fn to_julian(utc: DateTime<Utc>) -> f64 {
    let seconds_since_epoch: i64 = utc.timestamp();

    ((seconds_since_epoch as f64) / 86400.0) + 2440587.5
}

/// Convert Julian-Days to UTC time-point
///
/// Julian-Days is the number of days (incl. fraction) since January 1, 4713 BC, noon, without the
/// leap seconds.
///
/// # Arguments
///
/// * `jd` - Julian days since January 1, 4713 BC noon (12:00)
///
fn to_utc(jd: f64) -> DateTime<Utc> {
    let secs_since_epoch = (jd - 2440587.5) * 86400.0;
    let nanos = secs_since_epoch * (1000.0 * 1000.0 * 1000.0);
    Utc.timestamp_nanos(nanos as i64)
}

/// Projecting value into range [0,..,PI]
///
/// # Arguments
///
/// * `x` - radiant, not normalized to range [-PI..PI]
fn in_pi(x: f64) -> f64 {
    let n = (x / PI2) as i64;
    let result = x - (n as f64 * PI2);
    if result < 0.0 {
        result + PI2
    } else {
        result
    }
}

/// Returns the eccentricity of earth ellipse
///
/// # Arguments
///
/// * `t` - time according to standard equinox J2000.0
///
fn eps(t: f64) -> f64 {
    return RAD * (23.43929111 + ((-46.8150) * t - 0.00059 * t * t + 0.001813 * t * t * t) / 3600.0);
}

/// Calculates equation of time, returning the tuple (delta-ascension, declination)
///
/// # Arguments
///
/// * `t` - time according to standard equinox J2000.0
fn berechne_zeitgleichung<T: FloatOps>(t: f64) -> (f64, f64) {
    let mut ra_mittel = 18.71506921 + 2400.0513369 * t + (2.5862e-5 - 1.72e-9 * t) * t * t;

    let m = in_pi(PI2 * (0.993133 + 99.997361 * t));
    let l = in_pi(
        PI2 * (0.7859453
            + m / PI2
            + (6893.0 * T::sin(m) + 72.0 * T::sin(2.0 * m) + 6191.2 * t) / 1296.0e3),
    );
    let e = eps(t);
    let mut ra = T::atan(T::tan(l) * T::cos(e));

    if ra < 0.0 {
        ra += PI;
    }
    if l > PI {
        ra += PI;
    }

    ra = 24.0 * ra / PI2;

    let dk = T::asin(T::sin(e) * T::sin(l));

    // Damit 0<=RA_Mittel<24
    ra_mittel = 24.0 * in_pi(PI2 * ra_mittel / 24.0) / PI2;

    let mut d_ra = ra_mittel - ra;
    if d_ra < -12.0 {
        d_ra += 24.0;
    }
    if d_ra > 12.0 {
        d_ra -= 24.0;
    }

    d_ra = d_ra * 1.0027379;

    return (d_ra, dk);
}

/// Returning Sunrise and Sunset (or PolarNight/PolarDay) at geo-pos `lat/lon` at time `t` (UTC)
///
/// # Arguments
///
/// * `utc` - UTC time-point (DateTime<Utc>)
/// * `lat` - latitude in WGS84 system, ranging from -90.0 to 90.0.
/// * `lon` - longitude in WGS84 system, ranging from -180.0 to 180.0
///
/// North of ca 67.4 degree and south of ca -67.4 it may be permanent night or day. In this
/// case may be SunriseAndSet::PolarNight or SunriseAndSet::PolarDay.
///
/// In case latitude or longitude are not in valid ranges, the function will return Err(BadParam).
///
/// Algorithm ported to Rust from  [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
/// Its accuracy is in the range of a few minutes.
pub fn sunrise_and_set<F: FloatOps>(
    utc: DateTime<Utc>,
    lat: f64,
    lon: f64,
) -> Result<SunriseAndSet, SpaError> {
    if -90.0 > lat || 90.0 < lat || -180.0 > lon || 180.0 < lon {
        return Err(SpaError::BadParam);
    }

    let jd = to_julian(utc);
    let t = (jd - JD2000) / 36525.0;
    let h = -50.0 / 60.0 * RAD;
    let b = lat * RAD; // geographische Breite
    let geographische_laenge = lon;

    let (ra_d, dk) = berechne_zeitgleichung::<F>(t);

    let aux = (F::sin(h) - F::sin(b) * F::sin(dk)) / (F::cos(b) * F::cos(dk));
    if aux >= 1.0 {
        Result::Ok(SunriseAndSet::PolarNight)
    } else if aux <= -1.0 {
        Result::Ok(SunriseAndSet::PolarDay)
    } else {
        let zeitdifferenz = 12.0 * F::acos(aux) / PI;

        let aufgang_lokal = 12.0 - zeitdifferenz - ra_d;
        let untergang_lokal = 12.0 + zeitdifferenz - ra_d;
        let aufgang_welt = aufgang_lokal - geographische_laenge / 15.0;
        let untergang_welt = untergang_lokal - geographische_laenge / 15.0;
        let jd_start = F::trunc(jd); // discard fraction of day

        let aufgang_jd = (jd_start as f64) - 0.5 + (aufgang_welt / 24.0);
        let untergang_jd = (jd_start as f64) - 0.5 + (untergang_welt / 24.0);

        //	let untergang_utc = untergang_lokal - geographische_laenge /15.0;
        let sunriseset = SunriseAndSet::Daylight(to_utc(aufgang_jd), to_utc(untergang_jd));
        Result::Ok(sunriseset)
    }
}

/// Returning solar position (azimuth and zenith-angle) at time `t` and geo-pos `lat` and `lon`
///
/// # Arguments
///
/// * `utc` - UTC time-point (DateTime<Utc>)
/// * `lat` - latitude in WGS84 system, ranging from -90.0 to 90.0.
/// * `lon` - longitude in WGS84 system, ranging from -180.0 to 180.0
///
/// In case latitude or longitude are not in valid ranges, the function will return Err(BadParam).
///
/// Algorithm ported to Rust from [http://www.psa.es/sdg/sunpos.htm](https://web.archive.org/web/20220308165815/http://www.psa.es/sdg/sunpos.htm)
/// The algorithm is accurate to within 0.5 minutes of arc for the year 1999 and following.
pub fn solar_position<F: FloatOps>(utc: DateTime<Utc>, lat: f64, lon: f64) -> Result<SolarPos, SpaError> {
    if -90.0 > lat || 90.0 < lat || -180.0 > lon || 180.0 < lon {
        return Err(SpaError::BadParam);
    }

    let decimal_hours =
        (utc.hour() as f64) + ((utc.minute() as f64) + (utc.second() as f64) / 60.0) / 60.0;

    // Calculate difference in days between the current Julian Day
    // and JD 2451545.0, which is noon 1 January 2000 Universal Time
    let elapsed_julian_days = to_julian(utc) - JD2000;

    // Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    // ecliptic in radians but without limiting the angle to be less than 2*Pi
    // (i.e., the result may be greater than 2*Pi)
    let (ecliptic_longitude, ecliptic_obliquity) = {
        let omega = 2.1429 - (0.0010394594 * elapsed_julian_days);
        let mean_longitude = 4.8950630 + (0.017202791698 * elapsed_julian_days); // Radians
        let mean_anomaly = 6.2400600 + (0.0172019699 * elapsed_julian_days);
        let ecliptic_longitude = mean_longitude
            + 0.03341607 * F::sin(mean_anomaly)
            + 0.00034894 * F::sin(2.0 * mean_anomaly)
            - 0.0001134
            - 0.0000203 * F::sin(omega);
        let ecliptic_obliquity =
            0.4090928 - 6.2140e-9 * elapsed_julian_days + 0.0000396 * F::cos(omega);
        (ecliptic_longitude, ecliptic_obliquity)
    };

    // Calculate celestial coordinates ( right ascension and declination ) in radians
    // but without limiting the angle to be less than 2*Pi (i.e., the result may be
    // greater than 2*Pi)
    let (declination, right_ascension) = {
        let sin_ecliptic_longitude = F::sin(ecliptic_longitude);
        let dy = F::cos(ecliptic_obliquity) * sin_ecliptic_longitude;
        let dx = F::cos(ecliptic_longitude);
        let mut right_ascension = F::atan2(dy, dx);
        if right_ascension < 0.0 {
            right_ascension = right_ascension + PI2;
        }
        let declination = F::asin(F::sin(ecliptic_obliquity) * sin_ecliptic_longitude);
        (declination, right_ascension)
    };

    // Calculate local coordinates ( azimuth and zenith angle ) in degrees
    let (azimuth, zenith_angle) = {
        let greenwich_mean_sidereal_time =
            6.6974243242 + 0.0657098283 * elapsed_julian_days + decimal_hours;
        let local_mean_sidereal_time = (greenwich_mean_sidereal_time * 15.0 + lon) * RAD;
        let hour_angle = local_mean_sidereal_time - right_ascension;
        let latitude_in_radians = lat * RAD;
        let cos_latitude = F::cos(latitude_in_radians);
        let sin_latitude = F::sin(latitude_in_radians);
        let cos_hour_angle = F::cos(hour_angle);
        let mut zenith_angle = F::acos(
            cos_latitude * cos_hour_angle * F::cos(declination)
                + F::sin(declination) * sin_latitude,
        );
        let dy = -F::sin(hour_angle);
        let dx = F::tan(declination) * cos_latitude - sin_latitude * cos_hour_angle;
        let mut azimuth = F::atan2(dy, dx);
        if azimuth < 0.0 {
            azimuth = azimuth + PI2;
        }
        azimuth = azimuth / RAD;
        // Parallax Correction
        let parallax = (EARTH_MEAN_RADIUS / ASTRONOMICAL_UNIT) * F::sin(zenith_angle);
        zenith_angle = (zenith_angle + parallax) / RAD;
        (azimuth, zenith_angle)
    };

    let solpos = SolarPos {
        azimuth: azimuth,
        zenith_angle: zenith_angle,
    };

    Result::Ok(solpos)
}

#[cfg(test)]
mod tests {
    use chrono::{TimeZone, Timelike, Datelike, Utc};

    use super::FloatOps;
    use super::berechne_zeitgleichung;
    use super::solar_position;
    use super::sunrise_and_set;
    use super::eps;
    use super::in_pi;
    use super::to_julian;
    use super::to_utc;
    use super::SunriseAndSet;
    use super::JD2000;
    use super::PI2;
    use core::f64::consts::FRAC_PI_2;
    use core::f64::consts::PI;

    // for testing bind the `std` library
    extern crate std;

    const LAT_DEG: f64 = 48.1;
    const LON_DEG: f64 = 11.6;


    // FloatOps for the std environment
    pub struct StdFloatOps;

    // FloatOps for the std environment, mapping directly onto f64 operations
    impl FloatOps for StdFloatOps {
        fn sin(x: f64) -> f64 { x.sin() }
        fn cos(x: f64) -> f64 { x.cos() }
        fn tan(x: f64) -> f64 { x.tan() }
        fn asin(x: f64) -> f64 { x.asin() }
        fn acos(x: f64) -> f64 { x.acos() }
        fn atan(x: f64) -> f64 { x.atan() }
        fn atan2(y: f64, x: f64) -> f64 { y.atan2(x) }
        fn trunc(x: f64) -> f64 { x.trunc() }
    }

    #[test]
    fn test_pi() {
        assert!(FRAC_PI_2 < PI);
        assert!(in_pi(LAT_DEG) < PI2);
        assert!(in_pi(LON_DEG) < PI2);
    }

    #[test]
    fn test_datetime() {
        let unix_time_secs: i64 = 1128057420; // 2005-09-30T5:17:00Z
        let utc = Utc.timestamp_nanos(unix_time_secs * (1000 * 1000 * 1000));

        assert_eq!(utc.year(), 2005);
        assert_eq!(utc.month(), 9);
        assert_eq!(utc.day(), 30);
        assert_eq!(utc.hour(), 5);
        assert_eq!(utc.minute(), 17);
    }

    #[test]
    /// test-vector from [https://de.wikipedia.org/wiki/Sonnenstand](https://web.archive.org/web/20220712235611/https://de.wikipedia.org/wiki/Sonnenstand)
    fn test_julian_day() {

        //  6. August 2006 um 6 Uhr ut
        let dt = Utc.with_ymd_and_hms(2006, 8, 6, 6, 0, 0)
            .single().unwrap();

        let jd = to_julian(dt);
        assert_eq!(jd, 2453953.75);
        assert_eq!(jd - JD2000, 2408.75);

        assert_eq!(dt, to_utc(jd));
    }

    #[test]
    /// test-vector from [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
    fn test_zeitgleichung() {

        let exp_jd = 2453644.0;
        let exp_t = 0.057467488021902803;
        let exp_e = 0.40907976105657956;
        let exp_d_ra = 0.18539782794253773;
        let exp_dk = -0.05148602985190724;

        let dt = Utc.with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
            .single().unwrap();

        let jd = to_julian(dt);
        assert_eq!(exp_jd, jd);

        let t = (jd - JD2000) / 36525.0;
        assert_eq!(exp_t, t);

        assert_eq!(exp_e, eps(t));

        let (d_ra, dk) = berechne_zeitgleichung::<StdFloatOps>(t);

        assert_eq!(exp_d_ra, d_ra);
        assert_eq!(exp_dk, dk);
    }

    #[test]
    /// test-vector from [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
    fn test_sunrise_and_set() {

        let dt = Utc.with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
            .single().unwrap();

        // geo-pos near Frankfurt/Germany
        let lat = 50.0;
        let lon = 10.0;

        let sunriseandset = sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap();

        match sunriseandset {
            SunriseAndSet::PolarDay => assert!(false),
            SunriseAndSet::PolarNight => assert!(false),
            SunriseAndSet::Daylight(sunrise, sunset) => {
                assert_eq!(sunrise.year(), 2005);
                assert_eq!(sunrise.month(), 9);
                assert_eq!(sunrise.day(), 30);
                assert_eq!(sunrise.hour(), 5);
                assert_eq!(sunrise.minute(), 17);

                assert_eq!(sunset.year(), 2005);
                assert_eq!(sunset.month(), 9);
                assert_eq!(sunset.day(), 30);
                assert_eq!(sunset.hour(), 16);
                assert_eq!(sunset.minute(), 59);
            }
        }
    }

    #[test]
    /// test-vector from [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
    fn test_sunrise_and_set_polarday() {
        let dt = Utc.with_ymd_and_hms(2005, 6, 30, 12, 0, 0)
            .single().unwrap();

        // geo-pos in northern polar region
        let lat = 67.4;
        let lon = 10.0;

        let sunriseandset = sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap();

        match sunriseandset {
            SunriseAndSet::PolarDay => assert!(true),
            SunriseAndSet::PolarNight => assert!(false),
            _ => assert!(false),
        }
    }

    #[test]
    /// test-vector from [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
    fn test_sunrise_and_set_polarnight() {
        let dt = Utc.with_ymd_and_hms(2005, 6, 30, 12, 0, 0)
            .single().unwrap();

        // geo-pos in southern polar region
        let lat = -68.0;
        let lon = 10.0;

        let sunriseandset = sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap();

        match sunriseandset {
            SunriseAndSet::PolarDay => assert!(false),
            SunriseAndSet::PolarNight => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    /// test-vector from [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
    fn test_solar_position() {
        let exp_azimuth = 195.51003782406534;
        let exp_zenith_angle = 54.03653683638118;

        let dt = Utc.with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
            .single().unwrap();

        // geo-pos near Frankfurt/Germany
        let lat = 50.0;
        let lon = 10.0;

        let solpos = solar_position::<StdFloatOps>(dt, lat, lon).unwrap();

        assert_eq!(exp_azimuth, solpos.azimuth);
        assert_eq!(exp_zenith_angle, solpos.zenith_angle);
    }
}
