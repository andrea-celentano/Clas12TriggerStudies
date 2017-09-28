package clas12;

import org.jlab.geom.prim.Point3D;

public class ReconTOFHit extends RawTOFHit implements HitWithPositionEnergyTime  {

    private double energy, time; /* Recon objects */
    private Point3D p0; /* p0 in the SECTOR ref. frame */

    public ReconTOFHit(int _Sector, int _Panel, int _Paddle, int _Id, float _EnergyL, float _EnergyR, float _TimeL, float _TimeR, double energy, double time) {
        super(_Sector, _Panel, _Paddle, _Id, _EnergyL, _EnergyR, _TimeL, _TimeR);
        this.energy = energy;
        this.time = time;
        p0 = new Point3D();
    }

    public double getEnergy() {
        return energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public double getTime() {
        return time;
    }

    public void setTime(double time) {
        this.time = time;
    }

    public Point3D getP0() {
        return p0;
    }

    public void setP0(Point3D p0) {
        this.p0 = p0;
    }

    @Override
    public Point3D getPosition() {
        // TODO Auto-generated method stub
        return p0;
    }

}
