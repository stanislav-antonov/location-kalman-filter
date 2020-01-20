package kalmanfilter;

import com.sun.istack.internal.NotNull;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public abstract class FrameReferenceUtils {
    // WGS 84 datum constants
    // https://en.wikipedia.org/wiki/World_Geodetic_System#Defining_Parameters

    // Equatorial radius
    private final static double a = 6378137.0; // in meters
    // Flattening
    private final static double f = 1.0 / 298.257223563;
    // Polar semi-minor axis
    private final static double b = a * (1.0 - f); // in meters, 6356752.3142
    // Used in various places
    private final static double ba_squares_ratio = Math.pow(b, 2) / Math.pow(a, 2);
    // First eccentricity squared,
    private final static double e_square = 1.0 - ba_squares_ratio; // 6.69437999014 × 10^−3
    // private final static double e_square = (2.0 - f) * f;

    // http://mediatum.ub.tum.de/doc/1285843/393492.pdf
    @NotNull
    public static double[] enuToGeodetic(@NotNull double[] enuCoords,
                                         @NotNull Location locationReference) {
        final double[] ecefReferenceCoords = geodeticToEcef(locationReference);
        final double[] ecefCoords = enuToEcef(enuCoords, ecefReferenceCoords, locationReference);

        return ecefToGeodetic(ecefCoords);
    }

    public static double[] enuToEcef(@NotNull double[] enuCoords,
                                     @NotNull double[] ecefReferenceCoords,
                                     @NotNull Location locationReference) {
        final double x = enuCoords[0];
        final double y = enuCoords[1];
        final double z = enuCoords[2];

        final double x_r = ecefReferenceCoords[0];
        final double y_r = ecefReferenceCoords[1];
        final double z_r = ecefReferenceCoords[2];

        final double phi_r = Math.toRadians(locationReference.getLatitude()); // phi == lat
        final double lambda_r = Math.toRadians(locationReference.getLongitude()); // lambda == lon

        final double sin_phi_r = Math.sin(phi_r);
        final double cos_phi_r = Math.cos(phi_r);
        final double sin_lambda_r = Math.sin(lambda_r);
        final double cos_lambda_r = Math.cos(lambda_r);

        final RealMatrix M = new Array2DRowRealMatrix(new double[][] {
                {-sin_lambda_r, -sin_phi_r * cos_lambda_r, cos_phi_r * cos_lambda_r},
                {cos_lambda_r,  -sin_phi_r * sin_lambda_r, cos_phi_r * sin_lambda_r},
                {0,             cos_phi_r,                 sin_phi_r               }
        });

        final RealVector reference = new ArrayRealVector(new double[]
                {x_r, y_r, z_r}
        );

        final RealVector subject = new ArrayRealVector(new double[]
                {x, y, z}
        );

        final RealVector result = M.operate(subject).add(reference);

        return result.toArray();
    }

    // https://hal.archives-ouvertes.fr/hal-01704943/document
    @NotNull
    public static double[] ecefToGeodetic(@NotNull double[] ecefCoords) {
        final double x = ecefCoords[0];
        final double y = ecefCoords[1];
        final double z = ecefCoords[2];

        final double w_square = Math.pow(x, 2) + Math.pow(y, 2);
        final double l = e_square / 2.0;
        final double m = w_square / Math.pow(a, 2);
        final double n = Math.pow(z, 2) * (1.0 - e_square) / Math.pow(a, 2);
        final double p = (m + n - 4.0 * Math.pow(l, 2)) / 6.0;
        final double G = m * n * Math.pow(l, 2);
        final double H = 2 * Math.pow(p, 3) + G;
        final double C = Math.cbrt(H + G + 2.0 * Math.sqrt(H * G)) / Math.cbrt(2.0);
        final double i = - (2.0 * Math.pow(l, 2) + m + n) / 2.0;
        final double P = Math.pow(p, 2);
        final double beta = (i / 3.0) - C - (P / C);
        final double k = Math.pow(l, 2) * (Math.pow(l, 2) - m - n);
        final double t = Math.sqrt( Math.sqrt( Math.pow(beta, 2) - k ) - ((beta + i) / 2.0) )
                - Math.signum(m - n) * Math.sqrt(Math.abs(beta - i) / 2.0);

        final double F = Math.pow(t, 4) + 2.0 * i * Math.pow(t, 2) + 2.0 * l * (m - n) * t + k;
        final double dFdt = 4.0 * Math.pow(t, 3) + 4.0 * i * t + 2.0 * l * (m - n);
        final double delta_t = - F / dFdt;
        final double u = t + delta_t + l;
        final double v = t + delta_t - l;
        final double w = Math.sqrt(w_square);

        // in doubt, maybe should use atan() instead of atan2()
        final double phi = Math.atan2(z * u, w * v); // lat

        final double delta_w = w * (1.0 - 1.0 / u);
        final double delta_z = z * (1.0 - (1.0 - e_square) / v);

        final double h = Math.signum(u - 1.0) * Math.sqrt(Math.pow(delta_w, 2) + Math.pow(delta_z, 2));

        // here atan2() is correct for sure
        final double lambda = Math.atan2(y, x); // lon

        // {lat, lng, attitude}
        return new double[] {Math.toDegrees(phi), Math.toDegrees(lambda), h};
    }

    @NotNull
    public static double[] geodeticToEcef(@NotNull Location location) {
        return geodeticToEcef(new double[] {location.getLatitude(), location.getLongitude(), location.getAltitude()});
    }

    // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
    // https://github.com/substack/geodetic-to-ecef/blob/master/index.js
    @NotNull
    public static double[] geodeticToEcef(@NotNull double[] geodeticCoords) {
        double latitude = geodeticCoords[0];
        double longitude = geodeticCoords[1];
        double h = geodeticCoords[2];

        final double phi = Math.toRadians(latitude); // phi == lat
        final double lambda = Math.toRadians(longitude); // lambda == lon

        final double sin_phi = Math.sin(phi);
        final double cos_phi = Math.cos(phi);

        final double N_phi = a / Math.sqrt(1.0 - e_square * sin_phi * sin_phi);

        final double x = (N_phi + h) * cos_phi * Math.cos(lambda);
        final double y = (N_phi + h) * cos_phi * Math.sin(lambda);
        final double z = (ba_squares_ratio * N_phi + h) * sin_phi;

        return new double[] {x, y, z};
    }

    // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    @NotNull
    public static double[] ecefToEnu(@NotNull double[] ecefCoords,
                                     @NotNull double[] ecefReferenceCoords,
                                     @NotNull Location locationReference) {
        final double x_p = ecefCoords[0];
        final double y_p = ecefCoords[1];
        final double z_p = ecefCoords[2];

        final double x_r = ecefReferenceCoords[0];
        final double y_r = ecefReferenceCoords[1];
        final double z_r = ecefReferenceCoords[2];

        final double phi_r    = Math.toRadians(locationReference.getLatitude()); // phi == lat
        final double lambda_r = Math.toRadians(locationReference.getLongitude()); // lambda == lon

        final double sin_phi_r = Math.sin(phi_r);
        final double cos_phi_r = Math.cos(phi_r);
        final double sin_lambda_r = Math.sin(lambda_r);
        final double cos_lambda_r = Math.cos(lambda_r);

        final RealMatrix A = new Array2DRowRealMatrix(new double[][] {
                {-sin_lambda_r,             cos_lambda_r,              0        },
                {-sin_phi_r * cos_lambda_r, -sin_phi_r * sin_lambda_r, cos_phi_r},
                {cos_phi_r * cos_lambda_r, cos_phi_r * sin_lambda_r,   sin_phi_r}
        });

        final RealVector v = new ArrayRealVector(new double[]
                {x_p - x_r, y_p - y_r, z_p - z_r}
        );

        final RealVector result = A.operate(v);

        return result.toArray();
    }

    // http://www.epsg.org/Guidancenotes.aspx
    // https://www.iogp.org/wp-content/uploads/2019/09/373-07-02.pdf
    // https://hal.archives-ouvertes.fr/hal-01704943/document
    // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#Geodetic_to/from_ENU_coordinates
    @NotNull
    public static double[] geodeticToEnu(@NotNull Location location, @NotNull Location locationReference) {
        final double[] ecefCoords = geodeticToEcef(location);
        final double[] ecefReferenceCoords = geodeticToEcef(locationReference);

        return ecefToEnu(ecefCoords, ecefReferenceCoords, locationReference);
    }
}
