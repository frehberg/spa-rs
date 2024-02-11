// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0>.
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This example needs to be run with the `std` feature.

use chrono::{TimeZone, Utc};
use spa::{sunrise_and_set, StdFloatOps, SunriseAndSet};

fn main() {
    // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
    let dt = Utc
        .with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
        .single()
        .unwrap();

    // geo-pos near Frankfurt/Germany
    let lat = 50.0;
    let lon = 10.0;

    match sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap() {
        SunriseAndSet::Daylight(sunrise, sunset) => {
            println!("Sunrise and set: {} ----> {}", sunrise, sunset)
        }
        SunriseAndSet::PolarDay | SunriseAndSet::PolarNight => panic!(),
    }
}
