package clas12;

import org.jlab.geom.prim.Point3D;

public class CrossMatcherToReconFTOFHits extends CrossMatcherToFTOFHits{

   

    public CrossMatcherToReconFTOFHits(AnalysisClass ana, int layer) {
        super(ana, layer);
    }

    @Override
    protected double distanceProjCrossToHit(Point3D intersectCross, HitWithDCPositionEnergyTimeInfo hit) {
        Point3D pHit = hit.getPosition();
        double d = intersectCross.distance(pHit);
        return d;
    }


   

}
