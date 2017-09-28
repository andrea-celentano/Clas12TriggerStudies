package clas12;

import org.jlab.geom.prim.Point3D;



public class ECCluster implements HitWithPositionEnergyTime {

	public int sector, layer;
	public double energy, time;
	public Point3D p0;

	public ECCluster() {
		p0 = new Point3D();
	}

	public double getEnergy() {
		return energy;
	}

    @Override
    public Point3D getPosition() {
        // TODO Auto-generated method stub
        return p0;
    }

    @Override
    public double getTime() {
        // TODO Auto-generated method stub
        return time;
    }

}