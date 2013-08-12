(function (window){
    'use strict';

    /**
     * Class JTSK_Converter
     * @author Josef Zamrzla
     *
     * JTSK coordinates converter
     * This is a Javascript port of Pascal code that was originally published at
     * http://www.geospeleos.com/Mapovani/WGS84toSJTSK/WGS84toSJTSK.htm
     * Original author: Jakub Kerhat
     */

    var JTSK_Converter = function () {

        /**
         * Calculate distance between two points
         * @param x1
         * @param y1
         * @param x2
         * @param y2
         * @returns {*}
         */
        this.distPoints = function (x1, y1, x2, y2)
        {
            var dist = this.hypot(x1 - x2, y1 - y2);
            if (dist < this.EPS) {
                return 0;
            }

            return dist;
        };

        /**
         * Coordinates transformation
         * @param xs
         * @param ys
         * @param zs
         * @returns {Array}
         */
        this.transformCoords = function (xs, ys, zs)
        {
            // coeficients of transformation from WGS-84 to JTSK
            var dx = -570.69, dy = -85.69, dz = -462.84; // shift
            var wx = 4.99821/3600 * Math.PI / 180, wy = 1.58676/3600 * Math.PI / 180, wz = 5.2611/3600 * Math.PI / 180; // rotation
            var m  = -3.543e-6; // scale

            var xn = dx + (1 + m) * (+xs + wz * ys - wy * zs);
            var yn = dy + (1 + m) * (-wz * xs + ys + wx * zs);
            var zn = dz + (1 + m) * (+wy * xs - wx * ys + zs);

            return [xn, yn, zn];
        };

        // helper Math functions
        this.deg2rad = function (deg) {
            return (deg / 180) * Math.PI;
        };

        this.rad2deg = function (rad) {
            return rad / Math.PI * 180;
        };

        this.hypot = function (x, y) {
            return Math.sqrt(x * x + y * y) || 0;
        };

    };

    JTSK_Converter.prototype.EPS = 1e-4; // relative accuracy

    /**
     * Conversion from JTSK to WGS-84 (by iteration)
     * @param x
     * @param y
     * @returns {{lat: number, lon: number}}
     */
    JTSK_Converter.prototype.JTSKtoWGS84 = function (x, y)
    {
        if (!(x && y)) {
            return {'lat': 0, 'lon': 0};
        }

        var delta = 5;
        var latitude = 49;
        var longitude = 14;
        var steps = 0;

        var jtsk, v1, v2, v3, v4;

        do {
            jtsk = this.WGS84toJTSK(latitude - delta, longitude - delta);
            if (jtsk.x && jtsk.y) {
                v1 = this.distPoints(jtsk.x, jtsk.y, x, y);
            } else {
                v1 = 1e32;
            }

            jtsk = this.WGS84toJTSK(latitude - delta, longitude + delta);
            if (jtsk.x && jtsk.y) {
                v2 = this.distPoints(jtsk.x, jtsk.y, x, y);
            } else {
                v2 = 1e32;
            }

            jtsk = this.WGS84toJTSK(latitude + delta, longitude - delta);
            if (jtsk.x && jtsk.y) {
                v3 = this.distPoints(jtsk.x, jtsk.y, x, y);
            } else {
                v3 = 1e32;
            }

            jtsk = this.WGS84toJTSK(latitude + delta, longitude + delta);
            if (jtsk.x && jtsk.y) {
                v4 = this.distPoints(jtsk.x, jtsk.y, x, y);
            } else {
                v4 = 1e32;
            }

            if ((v1 <= v2) && (v1 <= v3) && (v1 <= v4)) {
                latitude = latitude - delta / 2;
                longitude = longitude - delta / 2;
            }

            if ((v2 <= v1) && (v2 <= v3) && (v2 <= v4)) {
                latitude = latitude - delta / 2;
                longitude = longitude + delta / 2;
            }

            if ((v3 <= v1) && (v3 <= v2) && (v3 <= v4)) {
                latitude = latitude + delta / 2;
                longitude = longitude - delta / 2;
            }

            if ((v4 <= v1) && (v4 <= v2) && (v4 <= v3)) {
                latitude = latitude + delta / 2;
                longitude = longitude + delta / 2;
            }

            delta *= 0.55;
            steps += 4;

        } while (!((delta < 0.00001) || (steps > 1000)));

        return {'lat': latitude, 'lon': longitude};
    };

    /**
     * Conversion from WGS-84 to JTSK
     * @param latitude
     * @param longitude
     * @returns {{x: number, y: number}}
     */
    JTSK_Converter.prototype.WGS84toJTSK = function (latitude, longitude)
    {
        if ((latitude < 40) || (latitude > 60) || (longitude < 5) || (longitude > 25)) {
            return {'x': 0, 'y': 0};
        } else {
            var lonlat = this.WGS84toBessel(latitude, longitude);
            return this.BesseltoJTSK(lonlat[0], lonlat[1]);
        }
    };

    /**
     * Conversion from ellipsoid WGS-84 to Bessel's ellipsoid
     * @param latitude
     * @param longitude
     * @param altitude
     * @returns {Array}
     */
    JTSK_Converter.prototype.WGS84toBessel = function (latitude, longitude, altitude)
    {
        var B = this.deg2rad(latitude);
        var L = this.deg2rad(longitude);
        var H = altitude || 0;

        var xyz1 = this.BLHToGeoCoords(B, L, H);
        var xyz2 = this.transformCoords(xyz1[0], xyz1[1], xyz1[2]);
        var BLH = this.geoCoordsToBLH(xyz2[0], xyz2[1], xyz2[2]);

        latitude = this.rad2deg(BLH[0]);
        longitude = this.rad2deg(BLH[1]);
        //Altitude = H;

        return [latitude, longitude];
    };

    /**
     * Conversion from Bessel's lat/lon to WGS-84
     * @param latitude
     * @param longitude
     * @returns {{x: number, y: number}}
     */
    JTSK_Converter.prototype.BesseltoJTSK = function (latitude, longitude)
    {
        var a     = 6377397.15508;
        var e     = 0.081696831215303;
        var n     = 0.97992470462083;
        var rho_0 = 12310230.12797036;
        var sinUQ = 0.863499969506341;
        var cosUQ = 0.504348889819882;
        var sinVQ = 0.420215144586493;
        var cosVQ = 0.907424504992097;
        var alfa  = 1.000597498371542;
        var k_2   = 1.00685001861538;

        var B = this.deg2rad(latitude);
        var L = this.deg2rad(longitude);

        var sinB = Math.sin(B);
        var t = (1 - e * sinB) / (1 + e * sinB);
        t = Math.pow(1 + sinB, 2) / (1 - Math.pow(sinB, 2)) * Math.exp(e * Math.log(t));
        t = k_2 * Math.exp(alfa * Math.log(t));

        var sinU  = (t - 1) / (t + 1);
        var cosU  = Math.sqrt(1 - sinU * sinU);
        var V     = alfa * L;
        var sinV  = Math.sin(V);
        var cosV  = Math.cos(V);
        var cosDV = cosVQ * cosV + sinVQ * sinV;
        var sinDV = sinVQ * cosV - cosVQ * sinV;
        var sinS  = sinUQ * sinU + cosUQ * cosU * cosDV;
        var cosS  = Math.sqrt(1 - sinS * sinS);
        var sinD  = sinDV * cosU / cosS;
        var cosD  = Math.sqrt(1 - sinD * sinD);

        var eps = n * Math.atan(sinD / cosD);
        var rho = rho_0 * Math.exp(-n * Math.log((1 + sinS) / cosS));

        return {'x': rho * Math.cos(eps), 'y': rho * Math.sin(eps)};
    };

    /**
     * Conversion from geodetic coordinates to Cartesian coordinates
     * @param B
     * @param L
     * @param H
     * @returns {Array}
     */
    JTSK_Converter.prototype.BLHToGeoCoords = function (B, L, H)
    {
        //  WGS-84 ellipsoid parameters
        var a   = 6378137.0;
        var f_1 = 298.257223563;
        var e2  = 1 - Math.pow(1 - 1 / f_1, 2);
        var rho = a / Math.sqrt(1 - e2 * Math.pow(Math.sin(B), 2));
        var x = (rho + H) * Math.cos(B) * Math.cos(L);
        var y = (rho + H) * Math.cos(B) * Math.sin(L);
        var z = ((1 - e2) * rho + H) * Math.sin(B);

        return [x, y, z];
    };

    /**
     * Conversion from Cartesian coordinates to geodetic coordinates
     * @param x
     * @param y
     * @param z
     * @returns {Array}
     */
    JTSK_Converter.prototype.geoCoordsToBLH = function (x, y, z)
    {
        // Bessel's ellipsoid parameters
        var a   = 6377397.15508;
        var f_1 = 299.152812853;
        var a_b = f_1 / (f_1-1);
        var p   = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2));
        var e2  = 1 - Math.pow(1 - 1 / f_1, 2);
        var th  = Math.atan(z * a_b / p);
        var st  = Math.sin(th);
        var ct  = Math.cos(th);
        var t   = (z + e2 * a_b * a * Math.pow(st, 3)) / (p - e2 * a * Math.pow(ct, 3));

        var B = Math.atan(t);
        var H = Math.sqrt(1 + t * t) * (p - a / Math.sqrt(1 + (1 - e2) * t * t));
        var L = 2 * Math.atan(y / (p + x));

        return [B, L, H];
    };

    window.JTSK_Converter = JTSK_Converter;

})(this);