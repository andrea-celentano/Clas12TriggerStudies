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

    private class Hit implements HitWithPositionEnergyTime {

        private Point3D p0;
        private double time, energy;
        private int detectorID; // 1-1A, 2-1B, 3-PCAL

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

        public Hit(double time, double energy, int detectorID) {
            super();
            this.time = time;
            this.energy = energy;
            this.detectorID = detectorID;
        }

        public Hit(ReconTOFHit hit, int i) {
            super();
            this.time = hit.getTime();
            this.energy = hit.getEnergy();
            this.detectorID = i;
            this.p0 = hit.getP0();
        }

        public Hit(ECCluster cluster, int i) {
            super();
            this.time = cluster.getTime();
            this.energy = cluster.getEnergy();
            this.detectorID = i;
            this.p0 = cluster.getPosition();
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

    private List<Hit> hits;
    Hit hitCurr;
    Hit hitNext;
    ListIterator<Hit> iteratorCurr;
    Integer posCurr;
    int nCoinc;
    
    private boolean debug;

    public HitsMatcherFTOF1PCal(AnalysisClass ana_, int sect_, int multiplicity_, double Tcoinc_, double Rcoinc_, double EminFTOF_, double EminPCAL_) {
        this.analysisClass = ana_;
        this.sector = sect_;
        this.Tcoinc = Tcoinc_;
        this.Rcoinc = Rcoinc_;
        this.multiplicity = multiplicity_;
        this.EminFTOF = EminFTOF_;
        this.EminPCAL = EminPCAL_;

        this.minPaddle1BAlone = 56;

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

        if (analysisClass.nevent==89583) debug=true;
        else debug=false;
        
        double Tmin, T0, T1, T2;
        T0 = 0;
        hits.clear();

        List<MatchedFTOF1PCALHits> theReturnList = new ArrayList<MatchedFTOF1PCALHits>();

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
                        theReturnList.add(new MatchedFTOF1PCALHits(hit1B.getTime(), this.sector, 1));
                    } else {
                        hits.add(new Hit(hit1B, 2));
                    }
                }
            }
        }
        if (!rejectPCAL) {
            for (ECCluster cluster : clustersPCAL) {
                if (cluster.getEnergy() > this.EminPCAL) {
                    hits.add(new Hit(cluster, 3));
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
        outerLoop: while (iteratorCurr.hasNext()) {

            posCurr++;
            hitCurr = iteratorCurr.next(); /* Get the base-hit (first in time window) */
            nCoinc = 1;
            T0 = hitCurr.getTime();

            if (debug) {
                System.out.println("HIT START BELOW");
                System.out.println(analysisClass.nevent + " " + (sector+1) + " -->" + T0 + " " + hitCurr.getDetectorID());
            }
            /* Easier case - multiplicity is 1. Don't need to do anything, just add the hit - unless it was required to HAVE another detector */
            if (multiplicity == 1) {
                theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, 1));
                continue;
            } else { /* For those cases when the multiplicity is > 1 */
                while (iteratorCurr.hasNext()) {
                    hitNext = iteratorCurr.next();
                    T1 = hitNext.getTime();

                    if (debug) System.out.println("NEXT TO COMPARE: " + T1 + " " + hitNext.getDetectorID());

                    if ((T1 - T0) > this.Tcoinc) {
                        iteratorCurr = hits.listIterator(posCurr);
                        if (nCoinc >= this.multiplicity) {
                            if (debug) System.out.println(" HIT END " + nCoinc);
                            theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, nCoinc));
                        }
                        continue outerLoop;
                    }
                    if (MatchHits(hitCurr, hitNext)) {
                        if (debug)  System.out.println("MATCHING HIT - REMOVE");
                        iteratorCurr.remove();
                        if (hitCurr.getDetectorID() != hitNext.getDetectorID()) nCoinc++;
                    }
                }
            }
          
            /* Repeat this because if that was the last hit, it won't call these */
            iteratorCurr = hits.listIterator(posCurr);            
            if (nCoinc >= this.multiplicity) {
                if (debug)  System.out.println(" HIT END " + nCoinc);
                theReturnList.add(new MatchedFTOF1PCALHits(T0, this.sector, nCoinc));
            }

        }

        return theReturnList;
    }

    private void SortHits() {
        Collections.sort(hits, new Comparator<HitWithPositionEnergyTime>() {
            @Override
            public int compare(HitWithPositionEnergyTime hit1, HitWithPositionEnergyTime hit2) {
                if (hit2.getTime() == hit1.getTime())
                    return 0;
                else if (hit1.getTime() < hit2.getTime())
                    return -1;
                else
                    return 1;
            }
        });
    }

    private boolean MatchHits(HitWithPositionEnergyTime hit1, HitWithPositionEnergyTime hit2) {
        boolean ret = false;
        if ((this.MatchGeo(hit1, hit2)) && (Math.abs(hit1.getTime() - hit2.getTime()) < this.Tcoinc)) ret = true;
        return ret;
    }

    private boolean MatchGeo(HitWithPositionEnergyTime hit1, HitWithPositionEnergyTime hit2) {
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
