# JTSK coordinates converter

Conversion example in PHP

```php
$converter = new JTSK\Converter();
$latlon = $converter->JTSKtoWGS84(1104335.13, 707849.38); // returns array ['lat', 'lon']
$jtsk = $converter->WGS84toJTSK(49.582556081943, 15.015852481984); // returns array ['x', 'y']
```

Conversion example in Javascript

```javascript
var converter = new JTSK_Converter();
var wgs = converter.JTSKtoWGS84(1104335.13, 707849.38); // returns object {'lat', 'lon'}
var jtsk = converter.WGS84toJTSK(49.582556081943, 15.015852481984); // returns object {'x', 'y'}
```