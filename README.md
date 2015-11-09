# JTSK coordinates converter


### Conversion example in Python

```python
converter = JtskConverter()
latlon = converter.jtsk_to_wgs84(1104335.13, 707849.38)  # returns dict {'lat', 'lon'}
jtsk = converter.wgs84_to_jtsk(49.582556081943, 15.015852481984)  # returns dict {'x', 'y'}
```

### Conversion example in Javascript

```javascript
var converter = new JTSK_Converter();
var wgs = converter.JTSKtoWGS84(1104335.13, 707849.38); // returns object {'lat', 'lon'}
var jtsk = converter.WGS84toJTSK(49.582556081943, 15.015852481984); // returns object {'x', 'y'}
```


### Conversion example in PHP

```php
$converter = new JTSK\Converter();
$latlon = $converter->JTSKtoWGS84(1104335.13, 707849.38); // returns array ['lat', 'lon']
$jtsk = $converter->WGS84toJTSK(49.582556081943, 15.015852481984); // returns array ['x', 'y']
```


## Live demo
Check a [live demo](http://vast-brushlands-3412.herokuapp.com/jtsk/) for examples of conversion accuracy
