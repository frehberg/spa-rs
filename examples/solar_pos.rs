// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0>.
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This example needs to be run with the `std` feature.

use chrono::{TimeZone, Utc};
use spa::{solar_position, SolarPos, StdFloatOps};

// main
fn main() {
    let dt = Utc
        .with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
        .single()
        .unwrap();

    // geo-pos near Frankfurt/Germany
    let lat = 50.0;
    let lon = 10.0;

    let solpos: SolarPos = solar_position::<StdFloatOps>(dt, lat, lon).unwrap();

    println!("Solar position: {}/{}", solpos.zenith_angle, solpos.azimuth);
}
