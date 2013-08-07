# JTSK coordinates converter

Conversion example

```php
$converter = new Converter();
$latlon = $converter->JTSKtoWGS84(1104335.13, 707849.38); // returns array ['lat', 'lon']
$jtsk = $converter->WGS84toJTSK(49.582556081943, 15.015852481984); // returns array ['x', 'y']
```