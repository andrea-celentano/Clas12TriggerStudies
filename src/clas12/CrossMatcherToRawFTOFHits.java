package clas12;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;

public class CrossMatcherToRawFTOFHits extends CrossMatcher {

    private int layerFTOF;

    /* Panel 1A FTOF - this is the back one */

    /* Panel 1B FTOF - this is the front one */
    private static double thetaAnglePlaneFTOF1B = 25;
    private static double thetaAngleMINFTOF1B = 3.667;
    private static double R2FTOF1B = 717.236;
    private static double counterWidthFTOF1B = 6;

    /* Panel 2 FTOT */
    private static double thetaAnglePlaneFTOF2 = 58.11;
    private static double thetaAngleMINFTOF2 = 34.698;
    private static double R2FTOF2 = 659.71;
    private static double counterWidthFTOF2 = 22;

    /* Following are key-parameters defining FTOF position */
    private double R2FTOF, thetaAnglePlaneFTOF, thetaAngleMINFTOF, counterWidthFTOF;

    /* Computed parameters */
    private double l0FTOF;
    private double hFTOF; // the distance between the
                          // (0,0,l0FTOF/cos(thetaAngleFTOF)) point and the
                          // start of FTOF along the FTOF-plane middle line
    private Plane3D planeFTOF;

    /* distances min values */
    private double minDistanceCSI, minDistanceY;

    public CrossMatcherToRawFTOFHits(AnalysisClass ana, int layer) {
        super(ana);
        layerFTOF = layer;
    }

    public double getMinDistanceCSI() {
        return minDistanceCSI;
    }

    public void setMinDistanceCSI(double minDistanceCSI) {
        this.minDistanceCSI = minDistanceCSI;
    }

    public double getMinDistanceY() {
        return minDistanceY;
    }

    public void setMinDistanceY(double minDistanceY) {
        this.minDistanceY = minDistanceY;
    }

    public double barLength(int barID) {
        /* See TOF-Geom note */
        double L = 0;
        switch (this.layerFTOF) {
        case 1:
            if (barID<=5) L = 15.85 * barID + 16.43;
            else L = 15.85 * barID +11.45;
            break;
        case 2:
            L = 6.40 * barID + 10.84;
            break;
        case 3:
            L = 13.73 * barID + 357.77;
            break;
        }
        return L;
    }

    public void setupGeo() {

        switch (this.layerFTOF) {
        case 1: /* Panel 1a-rear */
            System.out.println("not yet implemented");
            break;
        case 2: /* Panel 1b- front */
            thetaAnglePlaneFTOF = CrossMatcherToRawFTOFHits.thetaAnglePlaneFTOF1B;
            thetaAngleMINFTOF = CrossMatcherToRawFTOFHits.thetaAngleMINFTOF1B;
            R2FTOF = CrossMatcherToRawFTOFHits.R2FTOF1B;
            counterWidthFTOF = CrossMatcherToRawFTOFHits.counterWidthFTOF1B;
            break;
        case 3: /* Panel2 */

            /*
             * See my notes on 14/9/2017. Extracted from D. Carman description
             * of CLAS12 FTOF
             */
            thetaAnglePlaneFTOF = CrossMatcherToRawFTOFHits.thetaAnglePlaneFTOF2;
            thetaAngleMINFTOF = CrossMatcherToRawFTOFHits.thetaAngleMINFTOF2;
            R2FTOF = CrossMatcherToRawFTOFHits.R2FTOF2;
            counterWidthFTOF = CrossMatcherToRawFTOFHits.counterWidthFTOF2;
            break;
        }

        l0FTOF = this.R2FTOF * Math.cos(Math.toRadians(this.thetaAnglePlaneFTOF - this.thetaAngleMINFTOF));
        hFTOF = this.R2FTOF * Math.sin(Math.toRadians(this.thetaAngleMINFTOF)) / Math.cos(Math.toRadians(this.thetaAnglePlaneFTOF));

        Vector3D n = new Vector3D(Math.sin(Math.toRadians(this.thetaAnglePlaneFTOF)), 0., Math.cos(Math.toRadians(this.thetaAnglePlaneFTOF)));
        Point3D p = new Point3D(0, 0, this.l0FTOF / Math.cos(Math.toRadians(this.thetaAnglePlaneFTOF)));
        planeFTOF = new Plane3D(p, n);
    }

    /*
     * Given a cross (in the SECTOR ref. frame) and a list of hits in the FTOF
     * system, in that sector, return the minimum distance between the track -
     * projected to the FTOF plane - and hits. The Hit geometrical information
     * is basically contained in the "paddle" identifier. Return -1 if no hits
     * are present
     */
    public boolean matchCrossToFTOFHits(Cross cross, List<ReconTOFHit> hits) {

        boolean ret = false;
        if (hits.size() == 0) return false;

        /* Do the matching between this cross and the FTOFHits in this sector */
        int sector = cross.get_Sector();

        /*
         * Create a line passing by the cross point and with direction equal to
         * the cross direction
         */
        Point3D pCross = cross.get_Point(); /* passage point */
        Vector3D vCross = cross.get_Dir().toVector3D();/* Direction */
        Line3D rayCross = new Line3D(pCross, vCross);

        /* Determine the intersection between this line and the FTOF-plane */
        Point3D intersectCross = new Point3D();
        if ((planeFTOF.intersectionRay(rayCross, intersectCross) != 1)) {
            System.out.println("The given cross does not intersect TOF plane? event: " + analysisClass.nevent + " sector: " + cross.get_Sector() + " " + pCross.toString() + " " + vCross.toString());
            return false;
        }

        for (RawTOFHit hit : hits) {
            if (matchCrossToFTOFHit(intersectCross, hit)) {
                switch (layerFTOF){
                case 1:
                    break;
                case 2:
                    analysisClass.getHistogram2D("h2_FTOF1BEnergyMatched_LR").fill(hit.get_EnergyL(), hit.get_EnergyR());
                    break;
                case 3:
                    analysisClass.getHistogram2D("h2_FTOF2EnergyMatched_LR").fill(hit.get_EnergyL(), hit.get_EnergyR());
                    break;
                }
                
               
                ret = true;
            }
        }

        return ret;
    }

    /*
     * This is doing the merge between the Cross projected on the TOF plane and
     * the TOF hit
     */
    public boolean matchCrossToFTOFHit(Point3D crossProj, RawTOFHit hit) {
        /*
         * Get the cross coordinates in the TOF-plane system (see 14/9/2017
         * notes)
         */

        boolean ret = false;
        double y = crossProj.y();
        double csi = Math.sqrt(crossProj.x() * crossProj.x() + (crossProj.z() - l0FTOF / Math.cos(Math.toRadians(thetaAnglePlaneFTOF))) * (crossProj.z() - l0FTOF / Math.cos(Math.toRadians(thetaAnglePlaneFTOF))));

        /* The bar coordinates */
        int paddleN = hit.get_Paddle();
        double csiBar = hFTOF - counterWidthFTOF / 2 + paddleN * counterWidthFTOF;
        double csiBarMin = csiBar - counterWidthFTOF / 2;
        double csiBarMax = csiBar + counterWidthFTOF / 2;

        double yBar = 0;
        double yBarMin = yBar - barLength(paddleN) / 2;
        double yBarMax = yBar + barLength(paddleN) / 2;

        /* Now check */

        /* Case1: crossProj within counter area */
        if ((csi > csiBarMin) && (csi < csiBarMax) && (y > yBarMin) && (y < yBarMax)) {
            ret = true;
        }
        /* Case2: crossProj within csi but not y */
        else if ((csi > csiBarMin) && (csi < csiBarMax)) {
            if ((y > yBarMax) && ((y - yBarMax) < minDistanceY))
                ret = true;
            else if ((y < yBarMin) && ((yBarMin - y) < minDistanceY)) ret = true;
        }
        /* Case3: crossProj within y but not csi */
        else if ((y > yBarMin) && (y < yBarMax)) {
            if ((csi > csiBarMax) && ((csi - csiBarMax) < minDistanceCSI))
                ret = true;
            else if ((csi < csiBarMin) && ((csiBarMin - csi) < minDistanceCSI)) ret = true;

        }
        /* Case4: crossProj outside any dimension */
        else {
            double yBorder, csiBorder, distance;
            yBorder = 0;
            csiBorder = 0;
            if ((csi > csiBarMax) && (y > yBarMax)) {
                csiBorder = csiBarMax;
                yBorder = yBarMax;
            } else if ((csi > csiBarMax) && (y < yBarMin)) {
                csiBorder = csiBarMax;
                yBorder = yBarMin;
            } else if ((csi < csiBarMin) && (y < yBarMin)) {
                csiBorder = csiBarMin;
                yBorder = yBarMin;
            } else if ((csi < csiBarMin) && (y > yBarMax)) {
                csiBorder = csiBarMin;
                yBorder = yBarMax;
            }

            distance = ((csi - csiBorder) * (csi - csiBorder)) / (minDistanceCSI * minDistanceCSI);
            distance += ((y - yBorder) * (y - yBorder)) / (minDistanceY * minDistanceY);
            if (distance <= 1) ret = true;

        }

        return ret;

    }

}
