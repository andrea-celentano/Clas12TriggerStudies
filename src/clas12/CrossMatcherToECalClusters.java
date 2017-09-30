package clas12;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;

public class CrossMatcherToECalClusters extends CrossMatcher {

    /* EC */
    private static double thetaAngleEC = 25.;
    private static double l0PCAL = 697.7;
    private static double hPCAL = 385.1;
    private static double bPCAL = 394.;
    private static double dPCAL = 290.8;

    private static double l0ECIN;
    private static double l0ECOUT;

    private int idEC;
    private double l0EC, bEC, hEC, dEC;

    public CrossMatcherToECalClusters(AnalysisClass ana,int idEC) {
        super(ana);
        this.idEC = idEC;
    }

    protected void setupPlane() {
        switch (this.idEC) {
        case 0:
            this.l0EC = CrossMatcherToECalClusters.l0PCAL;
            this.bEC = CrossMatcherToECalClusters.bPCAL;
            this.hEC = CrossMatcherToECalClusters.hPCAL;
            this.dEC = CrossMatcherToECalClusters.dPCAL;
            break;
        case 1:
            this.l0EC = CrossMatcherToECalClusters.l0ECIN;
            break;
        case 2:
            this.l0EC = CrossMatcherToECalClusters.l0ECOUT;
            break;
        default:
            this.l0EC = CrossMatcherToECalClusters.l0PCAL;
            break;
        }

        Vector3D n = new Vector3D(Math.sin(Math.toRadians(CrossMatcherToECalClusters.thetaAngleEC)), 0., Math.cos(Math.toRadians(CrossMatcherToECalClusters.thetaAngleEC)));
        Point3D p = new Point3D(0, 0, this.l0EC / Math.cos(Math.toRadians(CrossMatcherToECalClusters.thetaAngleEC)));
        this.planeDetector = new Plane3D(p, n);
    }

    @Override
    protected double distanceProjCrossToHit(Point3D intersectCross, HitWithDCPositionEnergyTimeInfo hit) {
        Point3D pHit = hit.getPosition();
        double d = intersectCross.distance(pHit);
        return d;
    }

    
    /*
     * public boolean doesCrossMatchToClusters(Cross cross, List<ECCluster> clusters) {
     * double d = this.matchCrossToClusters(cross, clusters);
     * return this.distanceIsSmallerThanMin(d);
     * }
     */

}
