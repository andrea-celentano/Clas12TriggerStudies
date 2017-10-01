package clas12;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.jlab.geom.prim.Point3D;

/*This class handles the matching between FTOF1b and FTOF1a layers*/

public class HitsMatcherFTOF1PCal {

    private class Hit implements HitWithDCPositionEnergyTimeInfo {

        private Point3D p0;
        private double time, energy;
        private boolean isMatchedToR3Cross;
        private int detectorID; // 0-1B, 1-1A, 2-PCAL, 3-ECin, 4-ECout

        @Override
        public Point3D getPosition() {
            return p0;
        }

        public void setP0(Point3D p0) {
            this.p0 = p0;
        }

        @Override
        public double getTime() {
            return time;
        }

        public void setTime(double time) {
            this.time = time;
        }

        @Override
        public double getEnergy() {
            return energy;
        }

        public void setEnergy(double energy) {
            this.energy = energy;
        }

        public int getDetectorID() {
            return detectorID;
        }

        public void setDetectorID(int detectorID) {
            this.detectorID = detectorID;
        }

        public Hit(double time, double energy, int detectorID, boolean matched) {
            super();
            this.time = time;
            this.energy = energy;
            this.detectorID = detectorID;
            this.isMatchedToR3Cross = matched;
        }

        public Hit(ReconTOFHit hit, int i) {
            super();
            this.time = hit.getTime();
            this.energy = hit.getEnergy();
            this.detectorID = i;
            this.p0 = hit.getP0();
            this.isMatchedToR3Cross = hit.isMatchedToR3CrossProjection();
        }

        public Hit(ECCluster cluster, int i) {
            super();
            this.time = cluster.getTime();
            this.energy = cluster.getEnergy();
            this.detectorID = i;
            this.p0 = cluster.getPosition();
            this.isMatchedToR3Cross = cluster.isMatchedToR3CrossProjection();
        }

        @Override
        public boolean isMatchedToR3Segments() {
            return isMatchedToR3Cross;
        }

        @Override
        public boolean isMatchedToR3CrossProjection() {
            // TODO Auto-generated method stub
            return isMatchedToR3Cross;
        }

        @Override
        public double distanceR3CrossProjection() {
            // TODO Auto-generated method stub
            return 0;
        }

        @Override
        public void setIsMatchedToR3Segments(boolean val) {
            // TODO Auto-generated method stub

        }

        @Override
        public void setIsMatchedToR3CrossProjection(boolean val) {
            // TODO Auto-generated method stub

        }

        @Override
        public void setDistanceR3CrossProjection(double val) {
            // TODO Auto-generated method stub

        }

    }

    static class MatchedFTOF1PCALHits {
        double Time;
        int sector;
        int multiplicity;

        public double getTime() {
            return Time;
        }

        public void setTime(double time) {
            Time = time;
        }

        public int getSector() {
            return sector;
        }

        public void setSector(int sector) {
            this.sector = sector;
        }

        public int getMultiplicity() {
            return multiplicity;
        }

        public void setMultiplicity(int multiplicity) {
            this.multiplicity = multiplicity;
        }

        public MatchedFTOF1PCALHits(double time, int sector, int multiplicity) {
            super();
            Time = time;
            this.sector = sector;
            this.multiplicity = multiplicity;

        }
    }

    private AnalysisClass analysisClass;
    private int sector;
    private double Tcoinc;
    private double Rcoinc;
    private double EminFTOF;
    private double EminPCAL;

    private boolean rejectFTOF1A;
    private boolean rejectFTOF1B;
    private boolean rejectPCAL;

    private int multiplicity;
    private int minPaddle1BAlone;
    private int DC;

    private List<Hit> hits;
    Hit hitCurr;
    Hit hitNext;
    ListIterator<Hit> iteratorCurr;
    Integer posCurr;
    int nCoinc;
    int nCoincDC;

    private boolean layerHasHit[];
    private boolean layerHasHitDC[];

    private boolean debug;

    public HitsMatcherFTOF1PCal(AnalysisClass ana_, int sect_, int multiplicity_, double Tcoinc_, double Rcoinc_, double EminFTOF_, double EminPCAL_, int DC_) {
        this.analysisClass = ana_;
        this.sector = sect_;
        this.Tcoinc = Tcoinc_;
        this.Rcoinc = Rcoinc_;
        this.multiplicity = multiplicity_;
        this.EminFTOF = EminFTOF_;
        this.EminPCAL = EminPCAL_;
        this.minPaddle1BAlone = 58;
//        this.minPaddle1BAlone = 65;
        this.DC = DC_;

        this.layerHasHit = new boolean[5];
        this.layerHasHitDC = new boolean[5];

        this.hits = new ArrayList<Hit>();
    }

    public int getMultiplicity() {
        return multiplicity;
    }

    public boolean isRejectFTOF1A() {
        return rejectFTOF1A;
    }

    public void setRejectFTOF1A(boolean rejectFTOF1A) {
        this.rejectFTOF1A = rejectFTOF1A;
    }

    public boolean isRejectFTOF1B() {
        return rejectFTOF1B;
    }

    public void setRejectFTOF1B(boolean rejectFTOF1B) {
        this.rejectFTOF1B = rejectFTOF1B;
    }

    public boolean isRejectPCAL() {
        return rejectPCAL;
    }

    public void setRejectPCAL(boolean rejectPCAL) {
        this.rejectPCAL = rejectPCAL;
    }

    public void setReject(boolean rejectFTOF1A, boolean rejectFTOF1B, boolean rejectPCAL) {
        this.rejectFTOF1A = rejectFTOF1A;
        this.rejectFTOF1B = rejectFTOF1B;
        this.rejectPCAL = rejectPCAL;
    }

    public int getMinPaddle1BAlone() {
        return minPaddle1BAlone;
    }

    public void setMinPaddle1BAlone(int minPaddle1BAlone) {
        this.minPaddle1BAlone = minPaddle1BAlone;
    }

    /*
     * This performs the matching.
     */
    public List<MatchedFTOF1PCALHits> MatchHits(List<ReconTOFHit> hitsFTOF1A, List<ReconTOFHit> hitsFTOF1B, List<ECCluster> clustersPCAL) {

      
        
        double Tmin, T0, T1, T2;
        T0 = 0;
        List<MatchedFTOF1PCALHits> theReturnList = new ArrayList<MatchedFTOF1PCALHits>();

        hits.clear();

        if ((hitsFTOF1A.size() == 0) && (hitsFTOF1B.size() == 0) && (clustersPCAL.size() == 0)) return theReturnList;

        /* First, apply the energy cut /DC cut (TBD) separately to the 3 array of hits */
        if (!rejectFTOF1A) {
            for (ReconTOFHit hit1A : hitsFTOF1A) {
                if (hit1A.getEnergy() > this.EminFTOF) {
                    hits.add(new Hit(hit1A, 1));
                }
            }
        }
        /* Here do the 1-B work-around too for large-angle panels */
        if (!rejectFTOF1B) {
            for (ReconTOFHit hit1B : hitsFTOF1B) {
                if (hit1B.getEnergy() > this.EminFTOF) {
                    if ((hit1B.get_Paddle() > this.minPaddle1BAlone)) {
                        if ((DC <= 1) || (hit1B.isMatchedToR3CrossProjection())) {
                            theReturnList.add(new MatchedFTOF1PCALHits(hit1B.getTime(), this.sector, 1));
                        }
                    } else {
                        hits.add(new Hit(hit1B, 0));
                    }
                }
            }
        }
        if (!rejectPCAL) {
            for (ECCluster cluster : clustersPCAL) {
                if (cluster.getEnergy() > this.EminPCAL) {
                    hits.add(new Hit(cluster, 2));
                }
            }
        }
        if (hits.size() == 0) return theReturnList;

        /* Sort hits in time */
        SortHits();
        /* Check minimum time */
        Tmin = hits.get(0).getTime();
        /* Get iterator */
        iteratorCurr = hits.listIterator();

        /* Now start the "sliding window" */
        posCurr = 0;
        while (iteratorCurr.hasNext()) {

            for (int ii = 0; ii < layerHasHit.length; ii++) {
                layerHasHit[ii] = false;
                layerHasHitDC[ii] = false;
            }

            posCurr++;
            hitCurr = iteratorCurr.next(); /* Get the base-hit (first in time window) */
            T0 = hitCurr.getTime();
            layerHasHit[hitCurr.getDetectorID()] = true;
            layerHasHitDC[hitCurr.getDetectorID()] = hitCurr.isMatchedToR3CrossProjection();

            if (debug) {
                System.out.println("HIT START BELOW");
                System.out.println(analysisClass.nevent + " " + (sector + 1) + " -->" + T0 + " " + hitCurr.getDetectorID());
            }
            /* Easier case - multiplicity is 1. Don't need to do anything, just add the hit - unless it was required to HAVE another detector */
            if (multiplicity == 1) {
                if ((DC <= 1) || (hitCurr.isMatchedToR3CrossProjection())) {
                    theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, 1));
                }
                continue;
            } else { /* For those cases when the multiplicity is > 1 */
                while (iteratorCurr.hasNext()) {

                    hitNext = iteratorCurr.next();
                    T1 = hitNext.getTime();

                    if (debug) System.out.println("NEXT TO COMPARE: " + T1 + " " + hitNext.getDetectorID());

                    if ((T1 - T0) > this.Tcoinc) {
                        break;
                    }
                    if (MatchHits(hitCurr, hitNext)) {
                        if (debug) System.out.println("MATCHING HIT - REMOVE");
                        layerHasHit[hitNext.getDetectorID()] = true;
                        layerHasHitDC[hitNext.getDetectorID()] = hitNext.isMatchedToR3CrossProjection();
                        iteratorCurr.remove();
                    }
                }
            }

            /* We are here because: T1-T0 > this.Tcoinc OR because there are no more hits next */
            iteratorCurr = hits.listIterator(posCurr);
            nCoinc = 0;
            nCoincDC = 0;
            for (int ii = 0; ii < layerHasHit.length; ii++) {
                if (layerHasHit[ii]) nCoinc++;
                if (layerHasHitDC[ii]) nCoincDC++;
                if (debug) System.out.println("Layer: " + ii + " " + layerHasHit[ii] + " " + layerHasHitDC[ii]);
            }

            if (debug) {
                System.out.println("LAYERS COINC: " + nCoinc + " " + nCoincDC + " REQUIRED: " + this.multiplicity);
            }

            switch (DC) {
            case 0: /* No DC requirements */
            case 1: /* Only presence of DC-segments */
                if (nCoinc >= this.multiplicity) theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, nCoinc));
                break;
            case 2: /* At least 1 DC-matched */
                if ((nCoinc >= this.multiplicity) && (nCoincDC >= 1)) theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, nCoinc));
                break;
            case 3:/* ALL 1 DC-matched */
                if ((nCoinc >= this.multiplicity) && (nCoincDC >= this.multiplicity)) theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, nCoinc));
                break;
            }
        }

        return theReturnList;
    }

    private void SortHits() {
        Collections.sort(hits, new Comparator<HitWithDCPositionEnergyTimeInfo>() {
            @Override
            public int compare(HitWithDCPositionEnergyTimeInfo hit1, HitWithDCPositionEnergyTimeInfo hit2) {
                if (hit2.getTime() == hit1.getTime())
                    return 0;
                else if (hit1.getTime() < hit2.getTime())
                    return -1;
                else
                    return 1;
            }
        });
    }

    private boolean MatchHits(HitWithDCPositionEnergyTimeInfo hit1, HitWithDCPositionEnergyTimeInfo hit2) {
        boolean ret = false;
        if ((this.MatchGeo(hit1, hit2)) && (Math.abs(hit1.getTime() - hit2.getTime()) < this.Tcoinc)) ret = true;
        return ret;
    }

    private boolean MatchGeo(HitWithDCPositionEnergyTimeInfo hit1, HitWithDCPositionEnergyTimeInfo hit2) {
        boolean ret = false;
        double distance;
        
      
       
        /* These are in SECTOR ref. frame */
        Point3D p1 = new Point3D(hit1.getPosition());
        Point3D p2 = new Point3D(hit2.getPosition());

        /* Go to the TILTED ref. frame: in this way, I just consider X-Y and not Z (almost the same for the 3 detectors so close-by) */
        p1.rotateY(-1. * Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
        p2.rotateY(-1. * Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

        p1.setZ(0.);
        p2.setZ(0.);

        distance = p1.distance(p2);

        if (distance < this.Rcoinc) ret = true;

        if (debug) System.out.println("MatchGEo: " + p1.toString() + " " + p2.toString() + " " + distance + " " + ret);

        return ret;
    }

    private boolean MatchGeoRaw(ReconTOFHit hit1A, ReconTOFHit hit1B) {
        boolean ret = false;

        int id1A = hit1A.get_Paddle();
        int id1B = hit1B.get_Paddle();

        double coincMean = 3.2 + 2.44 * id1A;
        double coincSigma = 2.;

        int coincINF = (int) Math.floor(coincMean - 3 * coincSigma);
        int coincSUP = (int) Math.ceil(coincMean + 3 * coincSigma);

        if ((id1B >= coincINF) && (id1B <= coincSUP)) ret = true;

        return ret;
    }

}
