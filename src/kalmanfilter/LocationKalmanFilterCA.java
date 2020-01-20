package kalmanfilter;

import com.sun.istack.internal.NotNull;
import com.sun.istack.internal.Nullable;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

import static kalmanfilter.FrameReferenceUtils.enuToGeodetic;
import static kalmanfilter.FrameReferenceUtils.geodeticToEnu;

// Kalman Filter for Constant Acceleration Model (CA)
// See: https://github.com/balzer82/Kalman/blob/master/Kalman-Filter-CA-2.ipynb?create=1

public class LocationKalmanFilterCA extends KalmanFilter {

    // Filter step in milliseconds
    private final static long DELTA_T = 100;

    private final static double PROCESS_NOISE_SIGMA = 0.001;
    private final static double ACC_MEASUREMENT_SIGMA = 10.0;

    @Nullable
    private Location mLastLocation;

    @Nullable
    private Location mReferenceLocation;

    public LocationKalmanFilterCA() {
        super(initialErrorCovariance(), stateTransitionModel(), processNoiseModel(), null);
    }

    public long getStep() {
        return DELTA_T;
    }

    public void predict() {
        predict(null);
    }

    public void update(@NotNull Location location) {
        // we don't take Z coord into account
        location.setAltitude(0);

        mLastLocation = location;
        if (mReferenceLocation == null) {
            // Reference location for geodetic -> enu conversion
            mReferenceLocation = location;
        }

        double[] measurementInput = geodeticToEnu(location, mReferenceLocation);
        measurementInput = Arrays.copyOf(measurementInput, measurementInput.length - 1);

        final double[][] measurementModel = measurementModelForLocation();
        final double[][] measurementNoiseModel = measurementNoiseModelForLocation(location);
        update(measurementInput, measurementModel, measurementNoiseModel);
    }

    public void update(@NotNull double[] acceleration) {
        final double[] measurementInput = Arrays.copyOf(acceleration, acceleration.length);
        final double[][] measurementModel = measurementModelForAcceleration();
        final double[][] measurementNoiseModel = measurementNoiseModelForAcceleration();

        update(measurementInput, measurementModel, measurementNoiseModel);
    }

    @Nullable
    public Location locationEstimation() {
        Location location = null;
        if (mLastLocation != null && mReferenceLocation != null) {
            final double[] stateEstimation = getStateEstimation();
            final double[] estimatedEnuCoords = new double[] { stateEstimation[0], stateEstimation[1], 0 };
            final double[] estimatedGeodeticCoords = enuToGeodetic(estimatedEnuCoords, mReferenceLocation);

            location = new Location(mLastLocation);
            location.setLatitude(estimatedGeodeticCoords[0]);
            location.setLongitude(estimatedGeodeticCoords[1]);
        }

        return location;
    }

    @NotNull
    private static double[][] measurementModelForLocation() {
        return new double[][] {
                {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},

                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        };
    }

    @NotNull
    private static double[][] measurementNoiseModelForLocation(@NotNull Location location) {
        final double sigma2 = Math.pow(location.getAccuracy(), 2);
        return new double[][] {
                {sigma2, 0.0},

                {0.0, sigma2},
        };
    }

    @NotNull
    private static double[][] measurementModelForAcceleration() {
        return new double[][] {
                {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},

                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
        };
    }

    @NotNull
    private static double[][] measurementNoiseModelForAcceleration() {
        final double sigma2 = Math.pow(ACC_MEASUREMENT_SIGMA, 2);
        return new double[][] {
                {sigma2, 0.0},

                {0.0, sigma2},
        };
    }

    @NotNull
    private static double[][] initialErrorCovariance() {
        return new double[][] {
                {100.0, 0.0, 0.0, 0.0,   0.0, 0.0},

                {0.0, 100.0, 0.0, 0.0,   0.0, 0.0},

                {0.0, 0.0,   10.0, 0.0,  0.0, 0.0},

                {0.0, 0.0,   0.0, 10.0,  0.0, 0.0},

                {0.0, 0.0,   0.0, 0.0,   1.0, 0.0},

                {0.0, 0.0,   0.0, 0.0,   0.0, 1.0},
        };
    }

    @NotNull
    private static double[][] stateTransitionModel() {
        final double dt = DELTA_T / 1000.0;
        final double dt2 = 0.5 * Math.pow(dt, 2.0);
        return new double[][] {
                {1, 0, dt, 0, dt2, 0},

                {0, 1, 0, dt, 0, dt2},

                {0, 0, 1, 0,  dt, 0 },

                {0, 0, 0, 1,  0,  dt},

                {0, 0, 0, 0,  1,  0 },

                {0, 0, 0, 0,  0,  1 },
        };
    }

    @NotNull
    private static double[][] processNoiseModel() {
        final double dt = DELTA_T / 1000.0;
        final RealMatrix G = new Array2DRowRealMatrix(new double[][] {
                { 0.5 * Math.pow(dt, 2) },
                { 0.5 * Math.pow(dt, 2) },
                { dt  },
                { dt  },
                { 1.0 },
                { 1.0 },
        });

        final RealMatrix G_tr = G.copy().transpose();

        return G.multiply(G_tr).scalarMultiply(PROCESS_NOISE_SIGMA).getData();
    }
}
