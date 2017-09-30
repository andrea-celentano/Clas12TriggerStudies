package clas12;

import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;

public abstract class CrossMatcherToFTOFHits extends CrossMatcher {

    protected int layerFTOF;

    /* Panel 1A FTOF - this is the back one */
    protected static double thetaAnglePlaneFTOF1A = 25;
    protected static double thetaAngleMINFTOF1A = 5.453;
    protected static double R2FTOF1A = 726.689;
    protected static double counterWidthFTOF1A = 15.01;

    /* Panel 1B FTOF - this is the front one */
    protected static double thetaAnglePlaneFTOF1B = 25;
    protected static double thetaAngleMINFTOF1B = 3.667;
    protected static double R2FTOF1B = 717.236;
    protected static double counterWidthFTOF1B = 6;

    /* Panel 2 FTOT */
    protected static double thetaAnglePlaneFTOF2 = 58.11;
    protected static double thetaAngleMINFTOF2 = 34.698;
    protected static double R2FTOF2 = 659.71;
    protected static double counterWidthFTOF2 = 22;

    /* Following are key-parameters defining FTOF position */
    protected double R2FTOF, thetaAnglePlaneFTOF, thetaAngleMINFTOF, counterWidthFTOF;

    /* Computed parameters */
    protected double l0FTOF;
    protected double hFTOF; // the distance between the
                            // (0,0,l0FTOF/cos(thetaAngleFTOF)) point and the
                            // start of FTOF along the FTOF-plane middle line

    public CrossMatcherToFTOFHits(AnalysisClass ana, int layer) {    
        super(ana);
        this.layerFTOF=layer;
    }

    public double barLength(int barID) {
        /* See TOF-Geom note */
        double L = 0;
        switch (this.layerFTOF) {
        case 1:
            if (barID <= 5)
                L = 15.85 * barID + 16.43;
            else
                L = 15.85 * barID + 11.45;
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

    @Override
    protected void setupPlane() {
        switch (this.layerFTOF) {
        /*
         * See my notes on 14/9/2017. Extracted from D. Carman description
         * of CLAS12 FTOF
         */
        case 1: /* Panel 1a-rear */
            thetaAnglePlaneFTOF = CrossMatcherToRawFTOFHits.thetaAnglePlaneFTOF1A;
            thetaAngleMINFTOF = CrossMatcherToRawFTOFHits.thetaAngleMINFTOF1A;
            R2FTOF = CrossMatcherToRawFTOFHits.R2FTOF1A;
            counterWidthFTOF = CrossMatcherToRawFTOFHits.counterWidthFTOF1A;
            break;
        case 2: /* Panel 1b- front */
            thetaAnglePlaneFTOF = CrossMatcherToRawFTOFHits.thetaAnglePlaneFTOF1B;
            thetaAngleMINFTOF = CrossMatcherToRawFTOFHits.thetaAngleMINFTOF1B;
            R2FTOF = CrossMatcherToRawFTOFHits.R2FTOF1B;
            counterWidthFTOF = CrossMatcherToRawFTOFHits.counterWidthFTOF1B;
            break;
        case 3: /* Panel2 */
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
        this.planeDetector = new Plane3D(p, n);
        System.out.println("CrossMatcherToFTOFHits plane setup done: (layer:) " +this.layerFTOF+" "+this.planeDetector.toString());
    }

}
