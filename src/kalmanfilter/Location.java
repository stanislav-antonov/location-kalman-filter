package kalmanfilter;

import com.sun.istack.internal.NotNull;

public class Location {

    private double mAltitude;
    private double mLatitude;
    private double mLongitude;
    private double mAccuracy;

    public Location(@NotNull Location location) {
        mAccuracy = location.getAccuracy();
        mAltitude = location.getAltitude();
        mLongitude = location.getLongitude();
        mLatitude = location.getLatitude();
    }

    public double getLatitude() {
        return mLatitude;
    }

    public double getLongitude() {
        return mLongitude;
    }

    public double getAltitude() {
        return mAltitude;
    }

    public double getAccuracy() {
        return mAccuracy;
    }

    public void setLatitude(double value) {
        mLatitude = value;
    }

    public void setLongitude(double value) {
        mLongitude = value;
    }

    public void setAltitude(double value) {
        mAltitude = value;
    }

    public void setAccuracy(double value) {
        mAccuracy = value;
    }
}
