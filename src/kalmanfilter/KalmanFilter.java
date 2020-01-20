package kalmanfilter;

import com.sun.istack.internal.NotNull;
import com.sun.istack.internal.Nullable;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;

public class KalmanFilter {
    // Internal filter state
    //
    private class State {
        // Posteriori state estimate at given time
        @NotNull
        private RealVector x;

        // Posteriori error covariance matrix
        // (a measure of the estimated accuracy of the state estimate)
        @NotNull
        private RealMatrix P;

        State(@NotNull RealVector x, @NotNull RealMatrix P) {
            this.x = x;
            this.P = P;
        }
    }

    @NotNull
    private State mStateEstimation;

    // State-transition model
    @NotNull
    final private RealMatrix A;

    // State-transition model - transposed matrix
    @NotNull
    final private RealMatrix A_tr;

    // Covariance process noise
    @NotNull
    final private RealMatrix Q;

    // Control-input model
    @Nullable
    final private RealMatrix B;

    // Identity matrix
    @NotNull
    final private RealMatrix I;

    public KalmanFilter(@NotNull double[][] initialErrorCovariance,
                        @NotNull double[][] stateTransitionModel,
                        @NotNull double[][] processNoiseModel,
                        @Nullable double[][] controlInputModel
    ) {
        final RealMatrix P = new Array2DRowRealMatrix(initialErrorCovariance);
        if (!P.isSquare()) {
            throw new IllegalArgumentException("Error covariance matrix must be square");
        }

        final int size = P.getColumnDimension();
        final double[] x = new double[size];
        Arrays.fill(x,1.0D);

        mStateEstimation = new State(new ArrayRealVector(x), P);
        I = MatrixUtils.createRealIdentityMatrix(size);

        A = new Array2DRowRealMatrix(stateTransitionModel);
        A_tr = A.copy().transpose();
        Q = new Array2DRowRealMatrix(processNoiseModel);

        B = controlInputModel != null ? new Array2DRowRealMatrix(controlInputModel) : null;
    }

    /**
     * @param controlInput Control input vector
     */
    public void predict(@Nullable double[] controlInput) {
        RealVector u = null;
        if (controlInput != null) {
            if (B == null) {
                throw new IllegalStateException();
            }

            u = new ArrayRealVector(controlInput);
            if (u.getDimension() != B.getColumnDimension()) {
                throw new DimensionMismatchException(u.getDimension(), B.getColumnDimension());
            }
        }

        RealVector x = A.operate(mStateEstimation.x);
        if (u != null) {
            x = x.add(B.operate(u));
        }

        mStateEstimation.x = new ArrayRealVector(x);
        mStateEstimation.P = (A.multiply(mStateEstimation.P).multiply(A_tr)).add(Q);
    }

    /**
     * @param measurementInput Observed data vector that supposed to be filtered
     * @return
     */
    public void update(@NotNull double[] measurementInput,
                       @NotNull double[][] measurementModel,
                       @NotNull double[][] measurementNoiseModel) {
        // Measurement model
        final RealMatrix H = new Array2DRowRealMatrix(measurementModel);
        final RealMatrix H_tr = H.copy().transpose();

        // Covariance measurement noise
        final RealMatrix R = new Array2DRowRealMatrix(measurementNoiseModel);

        // Innovation (or pre-fit residual) covariance
        final RealMatrix S = (H.multiply(mStateEstimation.P).multiply(H_tr)).add(R);

        // Optimal Kalman gain
        final RealMatrix K = (mStateEstimation.P.multiply(H_tr)).multiply(MatrixUtils.inverse(S));

        // Measurement vector
        final RealVector z = new ArrayRealVector(measurementInput);

        // Innovation or measurement
        final RealVector y = z.subtract(H.operate(mStateEstimation.x));

        // Updated (a posteriori) state estimate
        mStateEstimation.x = mStateEstimation.x.add(K.operate(y));

        // Updated (a posteriori) estimate covariance
        mStateEstimation.P = (I.subtract(K.multiply(H))).multiply(mStateEstimation.P);
    }

    @NotNull
    public double[] getStateEstimation() {
        return mStateEstimation.x.toArray();
    }
}
