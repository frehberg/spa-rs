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
The SPA library supports both `std` and `no_std` target builds via the `FloatOps` trait:

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

- On `std` target builds (default), the built-in implementation `StdFloatOps` is provided
for your convenience, gated behind the `std` feature.
```toml
[dependencies]
spa = "^0.5"
```
- On `no_std` target builds, you need to provide your own implementation for
floating point operations, for example using [libm](https://docs.rs/libm/0.2.7/libm/).
The `default-features = false` option must be specified in a dependency declaration.
```toml
[dependencies]
spa = { version = "^0.5", default-features = false }
```

## Migrating from SPA v0.3 to SPA v0.5 and later

The version SPA 0.4 introduced the `std`/`no_std` abstraction `FloatOps` and functions have been renamed following Rust naming conventions.

| SPA v0.3              | SPA v0.5                       |
|-----------------------|--------------------------------|
|  `calc_sunrise_and_set` | `sunrise_and_set::<FloatOps>` |
|  `calc_solar_position`  | `solar_position::<FloatOps>` |

The version 0.5 defined the `std` implemenation `StdFloatOps` for `FloatOps`, and introduced dependency to crate `thiserror` for `std` targets.

Now, migrating `std` code from SPA v0.3 to SPA v0.5 is done, performing a renaming simply:   

Sample code using SPA v0.3
``` 
...
use spa::{calc_sunrise_and_set, calc_solar_position, SolarPos, SunriseAndSet, SpaError};

fn position(dt: DateTime<Utc>, lat: f64, lon: f64) -> Result<SolarPos, SpaError>  {
    calc_solar_position(dt, lat, lon)
}

fn rise_and_set(dt: DateTime<Utc>, lat: f64, lon: f64) -> Result<SunriseAndSet, SpaError>  {
    calc_sunrise_and_set(dt, lat, lon)
}
``` 

Sample code using SPA v0.5
``` 
...
use spa::{sunrise_and_set, solar_position, SolarPos, SunriseAndSet, SpaError, StdFloatOps};

fn position(dt: DateTime<Utc>, lat: f64, lon: f64) -> Result<SolarPos, SpaError>  {
    solar_position::<StdFloatOps>(dt, lat, lon)
}

fn rise_and_set(dt: DateTime<Utc>, lat: f64, lon: f64) -> Result<SunriseAndSet, SpaError>  {
    sunrise_and_set::<StdFloatOps>(dt, lat, lon)
}
``` 
