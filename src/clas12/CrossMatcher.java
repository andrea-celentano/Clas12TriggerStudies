package clas12;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;

public abstract class CrossMatcher {

    protected AnalysisClass analysisClass;
    protected Plane3D planeDetector;
    protected double minDistance;

    protected abstract void setupPlane();

    protected abstract double distanceProjCrossToHit(Point3D intersectCross, HitWithDCPositionEnergyTimeInfo hit);

    private String h1_matchedHitsE;

    private String h1_closerDistance;

    public CrossMatcher(AnalysisClass ana) {
        this.analysisClass = ana;
        h1_matchedHitsE = "";
        h1_closerDistance = "";

    }

    public void setupGeo() {
        this.setupPlane();
    }

    public double getMinDistance() {
        return minDistance;
    }

    public void setMinDistance(double minDistance) {
        this.minDistance = minDistance;
    }

    public void setAnalysisClass(AnalysisClass ana) {
        this.analysisClass = ana;
    }

    public String getH1_matchedHitsE() {
        return h1_matchedHitsE;
    }

    public void setH1_matchedHitsE(String h1_matchedHitsE) {
        this.h1_matchedHitsE = h1_matchedHitsE;
    }

    public String getH1_closerDistance() {
        return h1_closerDistance;
    }

    public void setH1_closerDistance(String h1_closerDistance) {
        this.h1_closerDistance = h1_closerDistance;
    }

    /*
     * Given a cross (in the SECTOR ref. frame) and a list of hits in that
     * sector, in the SECTOR ref. frame, return the minimum distance between the
     * track - projected to the EC plane - and ALL hits. Return -1 if no matching hits are
     * present.
     * Also, set the flag for the HitWithDCPositionEnergyTimeInfo for all the hits that are within the threshold distance with this cross.
     * Also, set the flag for the cross itself
     */
    public double matchCrossToHits(MatchedCross cross, List<? extends HitWithDCPositionEnergyTimeInfo> hits) {

        double minDist = -1;
        if (hits.size() == 0) return minDist;
        int sector = cross.get_Sector();
        double d;
        int imin = 0;
        ArrayList<Double> distances = new ArrayList<Double>();

        /*
         * Create a line passing by the cross point and with direction equal to the
         * cross direction
         */
        Point3D pCross = cross.get_Point(); /* passage point */
        Vector3D vCross = cross.get_Dir().toVector3D();/* Direction */
        Line3D rayCross = new Line3D(pCross, vCross);

        /* Determine the intersection between this line and the EC-plane */
        Point3D intersectCross = new Point3D();
        if ((planeDetector.intersectionRay(rayCross, intersectCross) != 1)) {
            System.out.println("HERE: The given cross does not intersect plane?" + pCross.toString() + " " + vCross.toString());
            return minDist;
        }

        /*
         * now, find the hits in this sector, and verify the matching, returning the distance.
         * If the return value is -1, then ignore this hit.
         * If one hit has been found, mark it. Also correct the minimum distance.
         */

        for (HitWithDCPositionEnergyTimeInfo hit : hits) {
            d = this.distanceProjCrossToHit(intersectCross, hit);

            if (d >= 0) {
                distances.add(new Double(d));
                if ((d < this.minDistance) || (d == 0.)) {
                    if ((!hit.isMatchedToR3CrossProjection() && (h1_matchedHitsE.length() > 0))) analysisClass.getHistogram1D(h1_matchedHitsE + "_" + sector).fill(hits.get(imin).getEnergy()); /* Just one time per hit */
                    hit.setIsMatchedToR3CrossProjection(true);
                    if (d < hit.distanceR3CrossProjection()) hit.setDistanceR3CrossProjection(d);
                }

            }
        }

        /* Find the minimum distance to associate it with the cross */
        minDist = 9999;
        for (Double dTmp : distances) {
            if (dTmp < minDist) {
                minDist = dTmp;
                imin = distances.indexOf(dTmp);
            }
        }

        if (h1_closerDistance.length() > 0) analysisClass.getHistogram1D(h1_closerDistance + "_" + sector).fill(minDist);

        return minDist;
    }

    public boolean distanceIsSmallerThanMin(double distance) {
        if ((distance < this.minDistance) && (distance > 0))
            return true;
        else
            return false;
    }

}
