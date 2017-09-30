package clas12;

import org.jlab.geom.prim.Point3D;

public class ECCluster implements HitWithDCPositionEnergyTimeInfo{

    public int sector, layer;
    public double energy, time;
    public Point3D p0;

    private boolean isMatchedToR3Segments;
    private boolean isMatchedToR3CrossProjection;
    private double distanceR3CrossProjection;

    public ECCluster() {
        p0 = new Point3D();
    }

    public double getEnergy() {
        return energy;
    }

    @Override
    public Point3D getPosition() {
        return p0;
    }

    @Override
    public double getTime() {
        return time;
    }

    @Override
    public boolean isMatchedToR3Segments() {
        return isMatchedToR3Segments;
    }

    @Override
    public boolean isMatchedToR3CrossProjection() {
        return isMatchedToR3CrossProjection;
    }

    @Override
    public double distanceR3CrossProjection() {
        return distanceR3CrossProjection;
    }

    @Override
    public void setIsMatchedToR3Segments(boolean val) {
        if (!val) {
            isMatchedToR3CrossProjection = false;
            isMatchedToR3Segments = false;
        } else {
            isMatchedToR3CrossProjection = true;

        }
    }

    @Override
    public void setIsMatchedToR3CrossProjection(boolean val) {
        if (val) {
            isMatchedToR3CrossProjection = true;
            isMatchedToR3Segments = true;
        } else {
            isMatchedToR3CrossProjection = false;
        }
    }

    @Override
    public void setDistanceR3CrossProjection(double val) {
        this.distanceR3CrossProjection=val;

    }

}