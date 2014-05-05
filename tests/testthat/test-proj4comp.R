context("proj4comp")

wgs84_a <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
wgs84_b <- "+init=epsg:4326"
test_that("wgs84 projection comparisons work", {
    expect_true(proj4comp(wgs84_a, wgs84_b))
})

utm16_full_epsg <- "+init=epsg:32616 +proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm16_full <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm16_simple <- "+proj=utm +zone=16 +ellps=WGS84 +units=m +no_defs"
utm16_epsg <- "+init=epsg:32616"
utm16_south <- "+init=epsg:32716"
test_that("utm projection comparisons work", {
    expect_true(proj4comp(utm16_full_epsg, utm16_full))
    expect_true(proj4comp(utm16_full_epsg, utm16_simple))
    expect_true(proj4comp(utm16_full_epsg, utm16_epsg))
    expect_true(proj4comp(utm16_full, utm16_simple))
    expect_true(proj4comp(utm16_full, utm16_epsg))
    expect_true(proj4comp(utm16_simple, utm16_epsg))
    expect_false(proj4comp("+init=epsg:32616", "+init=epsg:32617"))
    expect_true(proj4comp("+proj=utm +zone=16 +ellps=WGS84", "+proj=utm +zone=16 +ellps=WGS84"))
    expect_false(proj4comp("+proj=utm +zone=16 +ellps=WGS84", "+proj=utm +zone=16 +ellps=NAD83"))
    expect_false(proj4comp("+proj=utm +zone=16 +ellps=WGS84", "+proj=utm +zone=16 +ellps=wgs84"))
    expect_false(proj4comp("+proj=utm +zone=16 +ellps=WGS84", "+proj=utm +zone=16 +ellps=WGS85"))
    expect_error(proj4comp("+proj=utm +zone=16 +ellps=WGS84", "+proj=utm +zone=16"))
    expect_true(proj4comp("+proj=utm +zone=16 +ellps=WGS84", "+proj=utm +zone=16 +ellps=WGS84 +units=m"))
    expect_false(proj4comp("+proj=utm +zone=16 +ellps=WGS84 +units=in", "+proj=utm +zone=16 +ellps=WGS84 +units=m"))
    expect_false(proj4comp(utm16_full, utm16_south))
    expect_false(proj4comp(utm16_simple, utm16_south))
    expect_false(proj4comp(utm16_epsg, utm16_south))
    expect_error(proj4comp("+proj=utm +zone=16", "+proj=utm +zone=16"))
})
