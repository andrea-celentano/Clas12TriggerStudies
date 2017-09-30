package clas12;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.cross.Cross;

public class CrossMatcherToRawFTOFHits extends CrossMatcherToFTOFHits {

    /* distances min values */
    private double minDistanceCSI, minDistanceY;
    
    

    public CrossMatcherToRawFTOFHits(AnalysisClass ana, int layer) {
        super(ana, layer);
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

    /*
     * Get the cross coordinates in the TOF-plane system (see 14/9/2017
     * notes)
     */
    protected double distanceProjCrossToHit(Point3D crossProj, HitWithDCPositionEnergyTimeInfo hitBase) {
        
        RawTOFHit hit=(RawTOFHit)hitBase; /*is this working?*/

        double ret = -1;
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
        int flagCase = 0;

        /* Case1: crossProj within counter area */
        if ((csi > csiBarMin) && (csi < csiBarMax) && (y > yBarMin) && (y < yBarMax)) {
            flagCase = 1;
            ret = 0.;
        }
        /* Case2: crossProj within csi but not y */
        else if ((csi > csiBarMin) && (csi < csiBarMax)) {
            flagCase = 2;
            if ((y > yBarMax) && ((y - yBarMax) < minDistanceY))
                ret = (y-yBarMax);
            else if ((y < yBarMin) && ((yBarMin - y) < minDistanceY)) ret = (yBarMin - y);
        }
        /* Case3: crossProj within y but not csi */
        else if ((y > yBarMin) && (y < yBarMax)) {
            flagCase = 3;
            if ((csi > csiBarMax) && ((csi - csiBarMax) < minDistanceCSI))
                ret = (csi-csiBarMax);
            else if ((csi < csiBarMin) && ((csiBarMin - csi) < minDistanceCSI)) ret = (csiBarMin - csi);

        }
        /* Case4: crossProj outside any dimension */
        else {
            flagCase = 4;
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
            if (distance <= 1) ret = Math.sqrt(((csi - csiBorder) * (csi - csiBorder))+((y - yBorder) * (y - yBorder)));
        }
        return ret;
    }

}
