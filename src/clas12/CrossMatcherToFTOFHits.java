package clas12;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;

public class CrossMatcherToFTOFHits extends CrossMatcher {

	private int layerFTOF;

	/* Panel 2 FTOT */
	private static double thetaAngleFTOF2 = 58.11;
	private static double l0FTOF2 = 605.397; /*
												 * See my notes on 14/9/2017. Extracted from D. Carman description of CLAS12
												 * FTOF
												 */
	private double l0FTOF, thetaAngleFTOF;
	private Plane3D planeFTOF;

	public CrossMatcherToFTOFHits(AnalysisClass ana, int layer) {
		super(ana);
		layerFTOF = layer;
	}

	public void setupGeo() {

		switch (this.layerFTOF) {
		case 1:
		case 2:
			System.out.println("not yet implemented");
			break;
		case 3:
			l0FTOF = CrossMatcherToFTOFHits.l0FTOF2;
			thetaAngleFTOF = CrossMatcherToFTOFHits.thetaAngleFTOF2;

		}

		Vector3D n = new Vector3D(Math.sin(Math.toRadians(this.thetaAngleFTOF)), 0., Math.cos(Math.toRadians(this.thetaAngleFTOF)));
		Point3D p = new Point3D(0, 0, this.l0FTOF / Math.cos(Math.toRadians(this.thetaAngleFTOF)));
		planeFTOF = new Plane3D(p, n);
	}

	/*
	 * Given a cross (in the SECTOR ref. frame) and a list of hits in the FTOF
	 * system, in that sector, return the minimum distance between the track -
	 * projected to the FTOF plane - and hits. The Hit geometrical information is
	 * basically contained in the "paddle" identifier. Return -1 if no hits are
	 * present
	 */
	public double matchCrossToFTOFHits(Cross cross, List<SimpleTOFHit> hits) {

		double minDist = -1;
		double d;
		int imin = 0;
		ArrayList<Double> distances = new ArrayList<Double>();
		if (hits.size() == 0) return minDist;

		/* Do the matching between this cross and the FTOFHits in this sector */
		int sector = cross.get_Sector();

		/*
		 * Create a line passing by the cross point and with direction equal to the
		 * cross direction
		 */
		Point3D pCross = cross.get_Point(); /* passage point */
		Vector3D vCross = cross.get_Dir().toVector3D();/* Direction */
		Line3D rayCross = new Line3D(pCross, vCross);

		
		
		
		/* Determine the intersection between this line and the FTOF-plane */
		Point3D intersectCross = new Point3D();
		if ((planeFTOF.intersectionRay(rayCross, intersectCross) != 1)) {
			System.out.println("The given cross does not intersect TOF plane?" + pCross.toString() + " " + vCross.toString());
			return minDist;
		}
		
		/*
		 * now, find the clusters in this sector, and compute the distance
		 */

		for (SimpleTOFHit hit : hits) {
			matchCrossToFTOFHit(intersectCross,hit);
		}

		return 0;
	}
	
	/*This is doing the merge between the Cross projected on the TOF plane and the TOF hit*/
	public void matchCrossToFTOFHit(Point3D crossProj, SimpleTOFHit hit) {
	//	System.out.println(crossProj.toString());
	}
	
}
