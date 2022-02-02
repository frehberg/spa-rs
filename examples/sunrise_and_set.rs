// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0>.
// This file may not be copied, modified, or distributed
// except according to those terms.

use chrono::{TimeZone, Utc};
use spa::{calc_sunrise_and_set, SunriseAndSet};

fn main() {

    // test-vector from http://lexikon.astronomie.info/zeitgleichung/neu.html
    let dt = Utc.ymd(2005, 9, 30).and_hms(12, 0, 0);

    // geo-pos near Frankfurt/Germany
    let lat = 50.0;
    let lon = 10.0;

    match calc_sunrise_and_set(dt, lat, lon).unwrap() {
        SunriseAndSet::Daylight(sunrise, sunset) =>
            println!("Sunrise and set: {} ----> {}", sunrise, sunset),
        SunriseAndSet::PolarDay | SunriseAndSet::PolarNight => panic!(),
    }
}