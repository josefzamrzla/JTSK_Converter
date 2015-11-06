import math
import numpy


class JtskCoverter:

    EPS = 1.0e-4  # relative accuracy

    def __init__(self):
        pass

    # Conversion from JTSK to WGS-84 (by iteration)
    def jtsk_to_wgs84(self, x, y):
        if not x or not y:
            return {'lat': 0, 'lon': 0}

        delta = 5.0
        latitude = 49.0
        longitude = 14.0
        steps = 0

        while True:
            print steps, latitude, longitude
            jtsk = self.wgs84_to_jtsk(latitude - delta, longitude - delta)
            if jtsk['x'] and jtsk['y']:
                v1 = self.dist_points(jtsk['x'], jtsk['y'], x, y)
            else:
                v1 = 1.0e32

            jtsk = self.wgs84_to_jtsk(latitude - delta, longitude + delta)
            if jtsk['x'] and jtsk['y']:
                v2 = self.dist_points(jtsk['x'], jtsk['y'], x, y)
            else:
                v2 = 1.0e32

            jtsk = self.wgs84_to_jtsk(latitude + delta, longitude - delta)
            if jtsk['x'] and jtsk['y']:
                v3 = self.dist_points(jtsk['x'], jtsk['y'], x, y)
            else:
                v3 = 1.0e32

            jtsk = self.wgs84_to_jtsk(latitude + delta, longitude + delta)
            if jtsk['x'] and jtsk['y']:
                v4 = self.dist_points(jtsk['x'], jtsk['y'], x, y)
            else:
                v4 = 1.0e32

            if v1 <= v2 and v1 <= v3 and v1 <= v4:
                latitude -= (delta / 2.0)
                longitude -= (delta / 2.0)

            if v2 <= v1 and v2 <= v3 and v2 <= v4:
                latitude -= (delta / 2.0)
                longitude += (delta / 2.0)

            if v3 <= v1 and v3 <= v2 and v3 <= v4:
                latitude += (delta / 2.0)
                longitude -= (delta / 2.0)

            if v4 <= v1 and v4 <= v2 and v4 <= v3:
                latitude += (delta / 2.0)
                longitude += (delta / 2.0)

            delta *= 0.55
            steps += 4

            if delta < 0.00001 or steps > 1000:
                break

        return {'lat': latitude, 'lon': longitude}

    # Conversion from WGS-84 to JTSK
    def wgs84_to_jtsk(self, latitude, longitude):
        if latitude < 40 or latitude > 60 or longitude < 5 or longitude > 25:
            return {'x': 0, 'y': 0}
        else:
            (latitude, longitude) = self.wgs84_to_bessel(latitude, longitude)
            return self.bessel_to_jtsk(latitude, longitude)

    # Conversion from ellipsoid WGS-84 to Bessel's ellipsoid
    def wgs84_to_bessel(self, latitude, longitude, altitude=0.0):
        b = numpy.deg2rad(latitude)
        l = numpy.deg2rad(longitude)
        h = altitude

        (x1, y1, z1) = self.blht_to_geo_coords(b, l, h)
        (x2, y2, z2) = self.transform_coords(x1, y1, z1)
        (b, l, h) = self.geo_coords_to_blh(x2, y2, z2)

        latitude = numpy.rad2deg(b)
        longitude = numpy.rad2deg(l)
        # Altitude = h

        return [latitude, longitude]

    # Conversion from Bessel's lat/lon to WGS-84
    @staticmethod
    def bessel_to_jtsk(latitude, longitude):
        # a = 6377397.15508
        e = 0.081696831215303
        n = 0.97992470462083
        rho_0 = 12310230.12797036
        sin_uq = 0.863499969506341
        cos_uq = 0.504348889819882
        sin_vq = 0.420215144586493
        cos_vq = 0.907424504992097
        alfa = 1.000597498371542
        k_2 = 1.00685001861538

        b = numpy.deg2rad(latitude)
        l = numpy.deg2rad(longitude)

        sin_b = math.sin(b)
        t = (1 - e * sin_b) / (1 + e * sin_b)
        t = math.pow(1 + sin_b, 2) / (1 - math.pow(sin_b, 2)) * math.exp(e * math.log(t))
        t = k_2 * math.exp(alfa * math.log(t))

        sin_u = (t - 1) / (t + 1)
        cos_u = math.sqrt(1 - sin_u * sin_u)
        v = alfa * l
        sin_v = math.sin(v)
        cos_v = math.cos(v)
        cos_dv = cos_vq * cos_v + sin_vq * sin_v
        sin_dv = sin_vq * cos_v - cos_vq * sin_v
        sin_s = sin_uq * sin_u + cos_uq * cos_u * cos_dv
        cos_s = math.sqrt(1 - sin_s * sin_s)
        sin_d = sin_dv * cos_u / cos_s
        cos_d = math.sqrt(1 - sin_d * sin_d)

        eps = n * math.atan(sin_d / cos_d)
        rho = rho_0 * math.exp(-n * math.log((1 + sin_s) / cos_s))

        return {'x': rho * math.cos(eps), 'y': rho * math.sin(eps)}

    # Conversion from geodetic coordinates to Cartesian coordinates
    @staticmethod
    def blht_to_geo_coords(b, l, h):
        # WGS-84 ellipsoid parameters
        a = 6378137.0
        f_1 = 298.257223563
        e2 = 1 - math.pow(1 - 1 / f_1, 2)
        rho = a / math.sqrt(1 - e2 * math.pow(math.sin(b), 2))
        x = (rho + h) * math.cos(b) * math.cos(l)
        y = (rho + h) * math.cos(b) * math.sin(l)
        z = ((1 - e2) * rho + h) * math.sin(b)

        return [x, y, z]

    # Conversion from Cartesian coordinates to geodetic coordinates
    @staticmethod
    def geo_coords_to_blh(x, y, z):
        # Bessel's ellipsoid parameters
        a = 6377397.15508
        f_1 = 299.152812853
        a_b = f_1 / (f_1-1)
        p = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
        e2 = 1 - math.pow(1 - 1 / f_1, 2)
        th = math.atan(z * a_b / p)
        st = math.sin(th)
        ct = math.cos(th)
        t = (z + e2 * a_b * a * math.pow(st, 3)) / (p - e2 * a * math.pow(ct, 3))

        b = math.atan(t)
        h = math.sqrt(1 + t * t) * (p - a / math.sqrt(1 + (1 - e2) * t * t))
        l = 2.0 * math.atan(y / (p + x))

        return [b, l, h]

    # Distance between two points
    @staticmethod
    def dist_points(x1, y1, x2, y2):
        dist = math.hypot(x1 - x2, y1 - y2)
        if dist < JtskCoverter.EPS:
            return 0.0

        return dist

    # Coordinates transformation
    @staticmethod
    def transform_coords(xs, ys, zs):
        # coeficients of transformation from WGS-84 to JTSK
        dx = -570.69
        dy = -85.69
        dz = -462.84  # shift
        wx = 4.99821/3600 * math.pi / 180
        wy = 1.58676/3600 * math.pi / 180
        wz = 5.2611/3600 * math.pi / 180  # rotation
        m = -3.543e-6  # scale

        xn = dx + (1 + m) * (+xs + wz * ys - wy * zs)
        yn = dy + (1 + m) * (-wz * xs + ys + wx * zs)
        zn = dz + (1 + m) * (+wy * xs - wx * ys + zs)

        return [xn, yn, zn]
