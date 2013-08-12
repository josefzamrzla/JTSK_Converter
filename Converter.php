<?php
namespace JTSK;

/**
 * Class Converter
 * @package JTSK
 * @author Josef Zamrzla
 *
 * JTSK coordinates converter
 * This is a PHP port of Pascal code that was originally published at
 * http://www.geospeleos.com/Mapovani/WGS84toSJTSK/WGS84toSJTSK.htm
 * Original author: Jakub Kerhat
 */
class Converter {

    const EPS = 1e-4; // relative accuracy

    /**
     * Conversion from JTSK to WGS-84 (by iteration)
     *
     * @param $x
     * @param $y
     * @return array
     */
    public function JTSKtoWGS84($x, $y)
    {
        if (!($x && $y)) {
            return array('lat' => 0, 'lon' => 0);
        }

        $delta = 5;
        $latitude = 49;
        $longitude = 14;
        $steps = 0;

        do {
            $jtsk = $this->WGS84toJTSK($latitude - $delta, $longitude - $delta);
            if ($jtsk['x'] && $jtsk['y']) {
                $v1 = $this->distPoints($jtsk['x'], $jtsk['y'], $x, $y);
            } else {
                $v1 = 1e32;
            }

            $jtsk = $this->WGS84toJTSK($latitude - $delta, $longitude + $delta);
            if ($jtsk['x'] && $jtsk['y']) {
                $v2 = $this->distPoints($jtsk['x'], $jtsk['y'], $x, $y);
            } else {
                $v2 = 1e32;
            }

            $jtsk = $this->WGS84toJTSK($latitude + $delta, $longitude - $delta);
            if ($jtsk['x'] && $jtsk['y']) {
                $v3 = $this->distPoints($jtsk['x'], $jtsk['y'], $x, $y);
            } else {
                $v3 = 1e32;
            }

            $jtsk = $this->WGS84toJTSK($latitude + $delta, $longitude + $delta);
            if ($jtsk['x'] && $jtsk['y']) {
                $v4 = $this->distPoints($jtsk['x'], $jtsk['y'], $x, $y);
            } else {
                $v4 = 1e32;
            }

            if (($v1 <= $v2) && ($v1 <= $v3) && ($v1 <= $v4)) {
                $latitude = $latitude - $delta / 2;
                $longitude = $longitude - $delta / 2;
            }

            if (($v2 <= $v1) && ($v2 <= $v3) && ($v2 <= $v4)) {
                $latitude = $latitude - $delta / 2;
                $longitude = $longitude + $delta / 2;
            }

            if (($v3 <= $v1) && ($v3 <= $v2) && ($v3 <= $v4)) {
                $latitude = $latitude + $delta / 2;
                $longitude = $longitude - $delta / 2;
            }

            if (($v4 <= $v1) && ($v4 <= $v2) && ($v4 <= $v3)) {
                $latitude = $latitude + $delta / 2;
                $longitude = $longitude + $delta / 2;
            }

            $delta *= 0.55;
            $steps += 4;

        } while (!(($delta < 0.00001) || ($steps > 1000)));

        return array('lat' => $latitude, 'lon' => $longitude);

    }

    /**
     * Conversion from WGS-84 to JTSK
     *
     * @param $latitude
     * @param $longitude
     * @return array
     */
    public function WGS84toJTSK($latitude, $longitude)
    {
        if (($latitude < 40) || ($latitude > 60) || ($longitude < 5) || ($longitude > 25)) {
            return array('x' => 0, 'y' => 0);
        } else {
            list($latitude, $longitude) = $this->WGS84toBessel($latitude, $longitude);
            return $this->BesseltoJTSK($latitude, $longitude);
        }
    }

    /**
     * Conversion from ellipsoid WGS-84 to Bessel's ellipsoid
     *
     * @param $latitude
     * @param $longitude
     * @param int $altitude
     * @return array
     */
    public function WGS84toBessel($latitude, $longitude, $altitude = 0)
    {
        $B = deg2rad($latitude);
        $L = deg2rad($longitude);
        $H = $altitude;

        list($x1, $y1, $z1) = $this->BLHToGeoCoords($B, $L, $H);
        list($x2, $y2, $z2) = $this->transformCoords($x1, $y1, $z1);
        list($B, $L, $H) = $this->geoCoordsToBLH($x2, $y2, $z2);

        $latitude = rad2deg($B);
        $longitude = rad2deg($L);
        //$Altitude = $H;

        return array($latitude, $longitude);
    }

    /**
     * Conversion from Bessel's lat/lon to WGS-84
     *
     * @param $latitude
     * @param $longitude
     * @return array
     */
    public function BesseltoJTSK($latitude, $longitude)
    {
        $a     = 6377397.15508;
        $e     = 0.081696831215303;
        $n     = 0.97992470462083;
        $rho_0 = 12310230.12797036;
        $sinUQ = 0.863499969506341;
        $cosUQ = 0.504348889819882;
        $sinVQ = 0.420215144586493;
        $cosVQ = 0.907424504992097;
        $alfa  = 1.000597498371542;
        $k_2   = 1.00685001861538;

        $B = deg2rad($latitude);
        $L = deg2rad($longitude);

        $sinB = sin($B);
        $t = (1 - $e * $sinB) / (1 + $e * $sinB);
        $t = pow(1 + $sinB, 2) / (1 - pow($sinB, 2)) * exp($e * log($t));
        $t = $k_2 * exp($alfa * log($t));

        $sinU  = ($t - 1) / ($t + 1);
        $cosU  = sqrt(1 - $sinU * $sinU);
        $V     = $alfa * $L;
        $sinV  = sin($V);
        $cosV  = cos($V);
        $cosDV = $cosVQ * $cosV + $sinVQ * $sinV;
        $sinDV = $sinVQ * $cosV - $cosVQ * $sinV;
        $sinS  = $sinUQ * $sinU + $cosUQ * $cosU * $cosDV;
        $cosS  = sqrt(1 - $sinS * $sinS);
        $sinD  = $sinDV * $cosU / $cosS;
        $cosD  = sqrt(1 - $sinD * $sinD);

        $eps = $n * atan($sinD / $cosD);
        $rho = $rho_0 * exp(-$n * log((1 + $sinS) / $cosS));

        return array('x' => $rho * cos($eps), 'y' => $rho * sin($eps));
    }

    /**
     * Conversion from geodetic coordinates to Cartesian coordinates
     *
     * @param $B
     * @param $L
     * @param $H
     * @return array
     */
    public function BLHToGeoCoords($B, $L, $H)
    {
        // WGS-84 ellipsoid parameters
        $a   = 6378137.0;
        $f_1 = 298.257223563;
        $e2  = 1 - pow(1 - 1 / $f_1, 2);
        $rho = $a / sqrt(1 - $e2 * pow(sin($B), 2));
        $x = ($rho + $H) * cos($B) * cos($L);
        $y = ($rho + $H) * cos($B) * sin($L);
        $z = ((1 - $e2) * $rho + $H) * sin($B);

        return array($x, $y, $z);
    }

    /**
     * Conversion from Cartesian coordinates to geodetic coordinates
     *
     * @param $x
     * @param $y
     * @param $z
     * @return array
     */
    public function geoCoordsToBLH($x, $y, $z)
    {
        // Bessel's ellipsoid parameters
        $a   = 6377397.15508;
        $f_1 = 299.152812853;
        $a_b = $f_1 / ($f_1-1);
        $p   = sqrt(pow($x, 2) + pow($y, 2));
        $e2  = 1 - pow(1 - 1 / $f_1, 2);
        $th  = atan($z * $a_b / $p);
        $st  = sin($th);
        $ct  = cos($th);
        $t   = ($z + $e2 * $a_b * $a * pow($st, 3)) / ($p - $e2 * $a * pow($ct, 3));

        $B = atan($t);
        $H = sqrt(1 + $t * $t) * ($p - $a / sqrt(1 + (1 - $e2) * $t * $t));
        $L = 2 * atan($y / ($p + $x));

        return array($B, $L, $H);
    }

    /**
     * Distance between two points
     * @param $x1
     * @param $y1
     * @param $x2
     * @param $y2
     * @return float|int
     */
    private function distPoints($x1, $y1, $x2, $y2)
    {
        $dist = hypot($x1 - $x2, $y1 - $y2);
        if ($dist < self::EPS) {
            return 0;
        }

        return $dist;
    }

    /**
     * Coordinates transformation
     *
     * @param $xs
     * @param $ys
     * @param $zs
     * @return array
     */
    private function transformCoords($xs, $ys, $zs)
    {
        // coeficients of transformation from WGS-84 to JTSK
        $dx = -570.69; $dy = -85.69; $dz = -462.84; // shift
        $wx = 4.99821/3600 * pi() / 180; $wy = 1.58676/3600 * pi() / 180; $wz = 5.2611/3600 * pi() / 180; // rotation
        $m  = -3.543e-6; // scale

        $xn = $dx + (1 + $m) * (+$xs + $wz * $ys - $wy * $zs);
        $yn = $dy + (1 + $m) * (-$wz * $xs + $ys + $wx * $zs);
        $zn = $dz + (1 + $m) * (+$wy * $xs - $wx * $ys + $zs);

        return array($xn, $yn, $zn);
    }
}