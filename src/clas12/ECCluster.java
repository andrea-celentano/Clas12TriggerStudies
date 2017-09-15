package clas12;

import org.jlab.geom.prim.Point3D;



public class ECCluster {

	public int sector, layer;
	public double energy, time;
	public Point3D p0;

	public ECCluster() {
		p0 = new Point3D();
	}

	public double getEnergy() {
		return energy;
	}

}