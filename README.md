# Solar Position Algorithm (SPA)
The Solar Position Algorithm module (SPA) for Rust calculates the sunrise-sunset and azimuth and zenith-angle for
specific geo-position and time (UTC); for example for solar-panel-alignment or automotive.

## Calculating the Sunrise and Sunset 

The following function is calculating the sunrise-sunset at geo-pos `lat/lon` at time `t` (UTC). 

The algorithm has been ported to Rust from [http://lexikon.astronomie.info/zeitgleichung/neu.html](https://web.archive.org/web/20170812034800/http://lexikon.astronomie.info/zeitgleichung/neu.html)
Its accuracy is in the range of a few minutes.

#### Arguments

 * `utc` - UTC time-point (DateTime<Utc>)
 * `lat` - latitude in WGS84 system, ranging from -90.0 to 90.0.
 * `lon` - longitude in WGS84 system, ranging from -180.0 to 180.0

The function returns a result of type `SunriseAndSet`. 

```rust
pub enum SunriseAndSet {
    PolarNight,
    PolarDay,
    Daylight(DateTime<Utc>, DateTime<Utc>),
}
```

The polar night occurs in the northernmost and southernmost regions of the Earth when the night lasts
for more than 24 hours. This occurs only inside the polar circles. The opposite phenomenon, the
polar day, or midnight sun, occurs when the Sun stays above the horizon for more than 24 hours.
"Night" is understood as the center of the Sun being below a free horizon, represented by the variants
`SunriseAndSet::PolarNight` or `SunriseAndSet::PolarDay`.

Since the atmosphere bends the rays of the Sun, the polar day is longer than the polar night,
and the area that is affected by polar night is somewhat smaller than the area of midnight sun.
The polar circle is located at a latitude between these two areas, at the latitude of 
approximately 66.5 degrees. The function is approximating the atmospheric bend using a height
of 0,833333 degrees.

The variant `SunriseAndSet::Daylight(DateTime<Utc>, DateTime<Utc>)` represents time of sunrise and sunset.

In case latitude or longitude are not in valid ranges, the function will return `Result::Err(BadParam)`.


```rust
pub fn sunrise_and_set<F: FloatOps>(utc: DateTime<Utc>, lat: f64, lon: f64) -> Result<SunriseAndSet, SpaError> {..}
```

## Calculating the Solar Position

The following functions is calculating the solar position (azimuth and zenith-angle)
at time `t` and geo-pos `lat/lon`

The algorithm has been ported to Rust from [http://www.psa.es/sdg/sunpos.htm](https://web.archive.org/web/20220308165815/http://www.psa.es/sdg/sunpos.htm)
The algorithm is accurate to within 0.5 minutes of arc for the year 1999 and following.

#### Arguments

* `utc` - UTC time-point (DateTime<Utc>)
* `lat` - latitude in WGS84 system, ranging from -90.0 to 90.0.
* `lon` - longitude in WGS84 system, ranging from -180.0 to 180.0

The function returns a result of type `SolarPos`. 

```rust
pub struct SolarPos {
    // horizontal angle measured clockwise from a north base line or meridian
    pub azimuth: f64,
    // the angle between the zenith and the centre of the sun's disc
    pub zenith_angle: f64,
}
```

In case latitude or longitude are not in valid ranges, the function will return `Result::Err(BadParam)`.

```rust
pub fn solar_position<F: FloatOps>(utc: DateTime<Utc>, lat: f64, lon: f64) -> Result<SolarPos, SpaError> {..}
```

## Platform specific Float Operations 
The SPA library supports `std` and `no_std` target builds, but requires the implementation of the trigonometric 
operations etc. as defined by the trait `FloatOps`

```rust
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
```

The example below illustrates the usage for 
the `std` target build, but this example can be adapted using the no_std crate [libm](https://docs.rs/libm/0.2.7/libm/) implemenating
the floating point operaitons for targets without specific hardware support.

```rust

use chrono::{TimeZone, Utc};
use spa::{solar_position, sunrise_and_set, SolarPos, FloatOps, SunriseAndSet};


// FloatOps for the std environment
pub struct StdFloatOps;

// FloatOps for the std environment, mapping directly onto f64 operations
impl FloatOps for StdFloatOps {
    fn sin(x: f64) -> f64 { x.sin() }    fn cos(x: f64) -> f64 { x.cos() }     fn tan(x: f64) -> f64 { x.tan() }
    fn asin(x: f64) -> f64 { x.asin() }  fn acos(x: f64) -> f64 { x.acos() }   fn atan(x: f64) -> f64 { x.atan() }
    fn atan2(y: f64, x: f64) -> f64 { y.atan2(x) }
    fn trunc(x: f64) -> f64 { x.trunc() }
}

// main
fn main() {
    let dt = Utc.with_ymd_and_hms(2005, 9, 30, 12, 0, 0)
        .single().unwrap();

    // geo-pos near Frankfurt/Germany
    let lat = 50.0;
    let lon = 10.0;

    let _solpos: SolarPos = solar_position::<StdFloatOps>(dt, lat, lon).unwrap();
    // ...
    let _sunrise_set: SunriseAndSet =  sunrise_and_set::<StdFloatOps>(dt, lat, lon).unwrap();
    // ...
}
```
