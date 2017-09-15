package clas12;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;



public class CrossMatcherToECClusters extends CrossMatcher{

	/* EC */
	private static double thetaAngleEC = 25.;
	private static double l0PCAL = 697.7;
	private static double hPCAL = 385.1;
	private static double bPCAL = 394.;
	private static double dPCAL = 290.8;

	private static double l0ECIN;
	private static double l0ECOUT;

	private double minDistance;

	private int idEC;
	private double l0EC, bEC, hEC, dEC;
	private Plane3D planeEC;

	
	
	public CrossMatcherToECClusters(AnalysisClass ana,int idEC) {
		super(ana);
		this.idEC = idEC;
		this.analysisClass=ana;
	}

	public double getMinDistance() {
		return minDistance;
	}

	public void setMinDistance(double minDistance) {
		this.minDistance = minDistance;
	}
	
	public void setAnalysisClass(AnalysisClass ana){
		this.analysisClass = ana;
	}

	public void setupGeo() {
		switch (this.idEC) {
		case 0:
			this.l0EC = CrossMatcherToECClusters.l0PCAL;
			this.bEC = CrossMatcherToECClusters.bPCAL;
			this.hEC = CrossMatcherToECClusters.hPCAL;
			this.dEC = CrossMatcherToECClusters.dPCAL;
			break;
		case 1:
			this.l0EC = CrossMatcherToECClusters.l0ECIN;
			break;
		case 2:
			this.l0EC = CrossMatcherToECClusters.l0ECOUT;
			break;
		default:
			this.l0EC = CrossMatcherToECClusters.l0PCAL;
			break;
		}

		Vector3D n = new Vector3D(Math.sin(Math.toRadians(CrossMatcherToECClusters.thetaAngleEC)), 0.,
				Math.cos(Math.toRadians(CrossMatcherToECClusters.thetaAngleEC)));
		Point3D p = new Point3D(0, 0, this.l0EC / Math.cos(Math.toRadians(CrossMatcherToECClusters.thetaAngleEC)));
		planeEC = new Plane3D(p, n);
	}

	/*
	 * Given a cross (in the SECTOR ref. frame) and a list of clusters in that
	 * sector, in the SECTOR ref. frame, return the minimum distance between the
	 * track - projected to the EC plane - and hits. Return -1 if not clusters are
	 * present
	 */
	public double matchCrossToClusters(Cross cross, List<ECCluster> clusters) {

		double minDist = -1;
		double d;
		int imin = 0;
		ArrayList<Double> distances = new ArrayList<Double>();

		if (clusters.size() == 0) return minDist;

		/* Do the matching between this cross and the ECCLusters in this sector */
		int sector = cross.get_Sector();

		/*
		 * Create a line passing by the cross point and with direction equal to the
		 * cross direction
		 */
		Point3D pCross = cross.get_Point(); /* passage point */
		Vector3D vCross = cross.get_Dir().toVector3D();/* Direction */
		Line3D rayCross = new Line3D(pCross, vCross);

		/* Determine the intersection between this line and the EC-plane */
		Point3D intersectCross = new Point3D();
		if ((planeEC.intersectionRay(rayCross, intersectCross) != 1)) {
			System.out.println("The given cross does not intersect EC plane?" + pCross.toString() + " " + vCross.toString());
			return minDist;
		}
		
		/*
		 * now, find the clusters in this sector, and compute the distance
		 */

		for (ECCluster cluster : clusters) {
			Point3D pCluster = cluster.p0;
			d = pCluster.distance(intersectCross);
			distances.add(new Double(d));
		}

		/* Find the minimum distance */
		minDist = 9999;
		for (Double dTmp : distances) {
			if (dTmp < minDist) {
				minDist = dTmp;
				imin = distances.indexOf(dTmp);
			}
		}
		
		if (minDist < this.minDistance) {
			analysisClass.getHistogram("h1_closerClustersE_"+sector).fill(clusters.get(imin).energy);
		}
		analysisClass.getHistogram("h1_minCrossClusterDistance_"+sector).fill(minDist);

		return minDist;
	}

	public boolean distanceIsSmallerThanMin(double distance) {
		if (distance < this.minDistance)
			return true;
		else
			return false;
	}

	public boolean doesCrossMatchToClusters(Cross cross, List<ECCluster> clusters) {
		double d = this.matchCrossToClusters(cross, clusters);
		return this.distanceIsSmallerThanMin(d);
	}

}
