// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0>.
// This file may not be copied, modified, or distributed
// except according to those terms.

use chrono::{TimeZone, Utc};
use spa::{sunrise_and_set, FloatOps, SunriseAndSet};


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

fn main() {

    // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
    let dt = Utc.with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
        .single().unwrap();

    // geo-pos near Frankfurt/Germany
    let lat = 50.0;
    let lon = 10.0;

    match sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap() {
        SunriseAndSet::Daylight(sunrise, sunset) =>
            println!("Sunrise and set: {} ----> {}", sunrise, sunset),
        SunriseAndSet::PolarDay | SunriseAndSet::PolarNight => panic!(),
    }
}