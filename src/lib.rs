#![cfg_attr(feature = "benchmark", feature(test))]
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0>.
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Solar Position Algorithm module for Rust
//!
//! Collection of algorithms calculating sunrise/sunset and azimuth/zenith-angle.

extern crate chrono;

use chrono::prelude::Utc;
use chrono::DateTime;
use chrono::TimeZone;
use chrono::Timelike;
use std::error;
use std::f64::consts::PI;
use std::time::{SystemTime, UNIX_EPOCH};

const PI2: f64 = PI * 2.0;
const RAD: f64 = 0.017453292519943295769236907684886;
const EARTH_MEAN_RADIUS: f64 = 6371.01;
// In km
const ASTRONOMICAL_UNIT: f64 = 149597890.0;
// In km
const JD2000: f64 = 2451545.0;

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
impl std::fmt::Display for SpaError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match *self {
            SpaError::BadParam => write!(f, "Latitude or longitude are not within valid ranges."),
        }
    }
}

/// Implementing Error trait for SpaError
impl error::Error for SpaError {}

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
    let systime: SystemTime = utc.into();

    let seconds_since_epoch = systime.duration_since(UNIX_EPOCH).unwrap().as_secs();

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
    let secs = secs_since_epoch.trunc();
    let nanos = (secs_since_epoch - secs) * (1000.0 * 1000.0 * 1000.0);
    Utc.timestamp(secs as i64, nanos as u32)
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
        (result + PI2)
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
fn berechne_zeitgleichung(t: f64) -> (f64, f64) {
    let mut ra_mittel = 18.71506921 + 2400.0513369 * t + (2.5862e-5 - 1.72e-9 * t) * t * t;

    let m = in_pi(PI2 * (0.993133 + 99.997361 * t));
    let l = in_pi(
        PI2 * (0.7859453
            + m / PI2
            + (6893.0 * f64::sin(m) + 72.0 * f64::sin(2.0 * m) + 6191.2 * t) / 1296.0e3),
    );
    let e = eps(t);
    let mut ra = f64::atan(f64::tan(l) * f64::cos(e));

    if ra < 0.0 {
        ra += PI;
    }
    if l > PI {
        ra += PI;
    }

    ra = 24.0 * ra / PI2;

    let dk = f64::asin(f64::sin(e) * f64::sin(l));

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
/// Algorithm ported to Rust from  http://lexikon.astronomie.info/zeitgleichung/neu.html
/// Its accuracy is in the range of a few minutes.
pub fn calc_sunrise_and_set(
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

    let (ra_d, dk) = berechne_zeitgleichung(t);

    let aux = (f64::sin(h) - f64::sin(b) * f64::sin(dk)) / (f64::cos(b) * f64::cos(dk));
    if aux >= 1.0 {
        Result::Ok(SunriseAndSet::PolarNight)
    } else if aux <= -1.0 {
        Result::Ok(SunriseAndSet::PolarDay)
    } else {
        let zeitdifferenz = 12.0 * f64::acos(aux) / PI;

        let aufgang_lokal = 12.0 - zeitdifferenz - ra_d;
        let untergang_lokal = 12.0 + zeitdifferenz - ra_d;
        let aufgang_welt = aufgang_lokal - geographische_laenge / 15.0;
        let untergang_welt = untergang_lokal - geographische_laenge / 15.0;
        let jd_start = jd.trunc(); // discard fraction of day

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
/// Algorithm ported to Rust from http://www.psa.es/sdg/sunpos.htm
/// The algorithm is accurate to within 0.5 minutes of arc for the year 1999 and following.
pub fn calc_solar_position(utc: DateTime<Utc>, lat: f64, lon: f64) -> Result<SolarPos, SpaError> {
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
            + 0.03341607 * f64::sin(mean_anomaly)
            + 0.00034894 * f64::sin(2.0 * mean_anomaly)
            - 0.0001134
            - 0.0000203 * f64::sin(omega);
        let ecliptic_obliquity =
            0.4090928 - 6.2140e-9 * elapsed_julian_days + 0.0000396 * f64::cos(omega);
        (ecliptic_longitude, ecliptic_obliquity)
    };

    // Calculate celestial coordinates ( right ascension and declination ) in radians
    // but without limiting the angle to be less than 2*Pi (i.e., the result may be
    // greater than 2*Pi)
    let (declination, right_ascension) = {
        let sin_ecliptic_longitude = f64::sin(ecliptic_longitude);
        let dy = f64::cos(ecliptic_obliquity) * sin_ecliptic_longitude;
        let dx = f64::cos(ecliptic_longitude);
        let mut right_ascension = f64::atan2(dy, dx);
        if right_ascension < 0.0 {
            right_ascension = right_ascension + PI2;
        }
        let declination = f64::asin(f64::sin(ecliptic_obliquity) * sin_ecliptic_longitude);
        (declination, right_ascension)
    };

    // Calculate local coordinates ( azimuth and zenith angle ) in degrees
    let (azimuth, zenith_angle) = {
        let greenwich_mean_sidereal_time =
            6.6974243242 + 0.0657098283 * elapsed_julian_days + decimal_hours;
        let local_mean_sidereal_time = (greenwich_mean_sidereal_time * 15.0 + lon) * RAD;
        let hour_angle = local_mean_sidereal_time - right_ascension;
        let latitude_in_radians = lat * RAD;
        let cos_latitude = f64::cos(latitude_in_radians);
        let sin_latitude = f64::sin(latitude_in_radians);
        let cos_hour_angle = f64::cos(hour_angle);
        let mut zenith_angle = f64::acos(
            cos_latitude * cos_hour_angle * f64::cos(declination)
                + f64::sin(declination) * sin_latitude,
        );
        let dy = -f64::sin(hour_angle);
        let dx = f64::tan(declination) * cos_latitude - sin_latitude * cos_hour_angle;
        let mut azimuth = f64::atan2(dy, dx);
        if azimuth < 0.0 {
            azimuth = azimuth + PI2;
        }
        azimuth = azimuth / RAD;
        // Parallax Correction
        let parallax = (EARTH_MEAN_RADIUS / ASTRONOMICAL_UNIT) * f64::sin(zenith_angle);
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
    extern crate chrono;

    use self::chrono::prelude::Utc;
    use self::chrono::TimeZone;
    use super::berechne_zeitgleichung;
    use super::calc_solar_position;
    use super::calc_sunrise_and_set;
    use super::eps;
    use super::in_pi;
    use super::to_julian;
    use super::to_utc;
    use super::SunriseAndSet;
    use super::JD2000;
    use super::PI2;
    use std::f64::consts::FRAC_PI_2;
    use std::f64::consts::PI;
    use tests::chrono::Datelike;
    use tests::chrono::Timelike;

    const LAT_DEG: f64 = 48.1;
    const LON_DEG: f64 = 11.6;

    #[test]
    fn test_pi() {
        assert!(FRAC_PI_2 < PI);
        assert!(in_pi(LAT_DEG) < PI2);
        assert!(in_pi(LON_DEG) < PI2);
    }

    #[test]
    fn test_datetime() {
        let unix_time_secs = 1128057420; // 2005-09-30T5:17:00Z
        let utc = Utc.timestamp(unix_time_secs, 0);

        assert_eq!(utc.year(), 2005);
        assert_eq!(utc.month(), 9);
        assert_eq!(utc.day(), 30);
        assert_eq!(utc.hour(), 5);
        assert_eq!(utc.minute(), 17);
    }

    #[test]
    fn test_julian_day() {
        // test-vector from https://de.wikipedia.org/wiki/Sonnenstand

        //  6. August 2006 um 6 Uhr ut
        let dt = Utc.ymd(2006, 8, 6).and_hms(6, 0, 0);

        let jd = to_julian(dt);
        assert_eq!(jd, 2453953.75);
        assert_eq!(jd - JD2000, 2408.75);

        assert_eq!(dt, to_utc(jd));
    }

    #[test]
    fn test_zeitgleichung() {
        // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html

        let exp_jd = 2453644.0;
        let exp_t = 0.057467488021902803;
        let exp_e = 0.40907976105657956;
        let exp_d_ra = 0.18539782794253773;
        let exp_dk = -0.05148602985190724;

        let dt = Utc.ymd(2005, 9, 30).and_hms(12, 0, 0);

        let jd = to_julian(dt);
        assert_eq!(exp_jd, jd);

        let t = (jd - JD2000) / 36525.0;
        assert_eq!(exp_t, t);

        assert_eq!(exp_e, eps(t));

        let (d_ra, dk) = berechne_zeitgleichung(t);

        assert_eq!(exp_d_ra, d_ra);
        assert_eq!(exp_dk, dk);
    }

    #[test]
    fn test_calc_sunrise_and_set() {
        // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
        let dt = Utc.ymd(2005, 9, 30).and_hms(12, 0, 0);

        // geo-pos near Frankfurt/Germany
        let lat = 50.0;
        let lon = 10.0;

        let sunriseandset = calc_sunrise_and_set(dt, lat, lon).unwrap();

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
    fn test_calc_sunrise_and_set_polarday() {
        // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
        let dt = Utc.ymd(2005, 6, 30).and_hms(12, 0, 0);

        // geo-pos in northern polar region
        let lat = 67.4;
        let lon = 10.0;

        let sunriseandset = calc_sunrise_and_set(dt, lat, lon).unwrap();

        match sunriseandset {
            SunriseAndSet::PolarDay => assert!(true),
            SunriseAndSet::PolarNight => assert!(false),
            _ => assert!(false),
        }
    }

    #[test]
    fn test_calc_sunrise_and_set_polarnight() {
        // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
        let dt = Utc.ymd(2005, 6, 30).and_hms(12, 0, 0);

        // geo-pos in southern polar region
        let lat = -68.0;
        let lon = 10.0;

        let sunriseandset = calc_sunrise_and_set(dt, lat, lon).unwrap();

        match sunriseandset {
            SunriseAndSet::PolarDay => assert!(false),
            SunriseAndSet::PolarNight => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn test_calc_solar_position() {
        // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
        let exp_azimuth = 195.51003782406534;
        let exp_zenith_angle = 54.03653683638118;

        let dt = Utc.ymd(2005, 9, 30).and_hms(12, 0, 0);

        // geo-pos near Frankfurt/Germany
        let lat = 50.0;
        let lon = 10.0;

        let solpos = calc_solar_position(dt, lat, lon).unwrap();

        assert_eq!(exp_azimuth, solpos.azimuth);
        assert_eq!(exp_zenith_angle, solpos.zenith_angle);
    }
}

/// # Benchmarks
/// Execute the benchmarks entering command:
/// ```commandline
/// cargo bench --features benchmark
/// ```
#[cfg(all(feature = "benchmark", test))]
mod benchmark {
    extern crate chrono;
    extern crate test;

    use self::chrono::prelude::Utc;
    use self::chrono::TimeZone;
    use super::calc_solar_position;
    use super::calc_sunrise_and_set;
    use std::time::SystemTime;

    #[bench]
    fn bench_reference_read_systime(b: &mut test::Bencher) {
        b.iter(|| {
            // any random value, tp prevent the code is discarded by code-optimizer
            let now: u64 = SystemTime::now().elapsed().unwrap().as_secs();

            // return the value to prevent the code is marked as unused and discarded by optimizer
            now
        });
    }

    #[bench]
    fn bench_calc_solar_position(b: &mut test::Bencher) {
        b.iter(|| {
            // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
            let exp_azimuth = 195.51003782406534;
            let exp_zenith_angle = 54.03653683638118;

            let dt = Utc.ymd(2005, 9, 30).and_hms(12, 0, 0);

            // geo-pos near Frankfurt/Germany
            let lat = 50.0;
            let lon = 10.0;

            let solpos = calc_solar_position(dt, lat, lon).unwrap();

            assert_eq!(exp_azimuth, solpos.azimuth);
            assert_eq!(exp_zenith_angle, solpos.zenith_angle);

            // return the value to prevent the code is marked as unused and discarded by optimizer
            solpos
        });
    }

    #[bench]
    fn bench_calc_sunrise_and_set(b: &mut test::Bencher) {
        b.iter(|| {
            // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
            let dt = Utc.ymd(2005, 9, 30).and_hms(12, 0, 0);

            // geo-pos near Frankfurt/Germany
            let lat = 50.0;
            let lon = 10.0;

            let sunriseandset = calc_sunrise_and_set(dt, lat, lon).unwrap();

            // return the value to prevent the code is marked as unused and discarded by optimizer
            sunriseandset
        });
    }
}
