package clas12;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.clas.physics.RecEvent;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.ui.TCanvas;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.*;
import org.jlab.io.evio.*;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Vector3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Line3D;
import org.jlab.rec.dc.cross.Cross;

import java.util.List;
import java.util.ArrayList;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

@SuppressWarnings("unused")
public class AnalysisClass {

    private String inputFileName = "";
    private int nEvents = -1;

    private DataReader dataReader;

    /* Here goes common definitions */
    public static double phiAngle_CLAS12 = 60.;
    public static double thetaAngleDC_CLAS12 = 25.;
    public static int nSectors_CLAS12 = 6;

    /* Generated particles */
    int nGenParticles;
    List<MatchedParticle> genParticles;

    /* Reconstructed particles */
    int nRecParticles;
    List<MatchedParticle> recParticles;

    /* Tracks */
    List<TrackMatchedToGen> tracksDC;

    /* DC crosses */
    public static short crossIsAssociatedToTrack = 255;
    List<MatchedCross> crossesDC;

    /* DC segments */
    List<SimpleDCSegment> segmentsDC;
    boolean[] sectorHasDCsegments;

    /* ECAL Clusters */
    private int idEC = 0; // 0: PCAL, 1: EC-IN, 2: EC-OUT
    private double minClusterE_ECAL = 0.;
    List<ECCluster>[] clusters;

    /* TOF */
    private double minE_FTOF2 = 0.1;
    private double minE_FTOF1B = 0.1;
    private double minE_FTOF1A = 0.1;
    List<ReconTOFHit>[] hitsFTOF2;
    List<ReconTOFHit>[] hitsFTOF1B;
    List<ReconTOFHit>[] hitsFTOF1A;

    /* TOF "time history matching */
    double Tmin = 100;
    double Tcoinc = 100;

    HitsMatcherFTOF1PCal hitsMatcherFTOF1[];

    int nFTOF1PCAL_thisEvent[] = new int[AnalysisClass.nSectors_CLAS12];
    int nFTOF2_thisEvent[] = new int[AnalysisClass.nSectors_CLAS12];

    int nMatchingsTotalALL_thisEvent;
    int nMatchingsDiffSectorsALL_thisEvent; // this is between 0 and 6

    int nMatchingsTotalFTOF1_thisEvent;
    int nMatchingsDiffSectorsFTOF1_thisEvent; // this is between 0 and 6

    int nMatchingsTotalFTOF2_thisEvent;
    int nMatchingsDiffSectorsFTOF2_thisEvent; // this is between 0 and 6

    int nMatchings1FTOF_total = 0;
    int nMatchings1FTOF1_total = 0;
    int nMatchings1FTOF2_total = 0;
    int nMatchings1FTOF_DifferentSectors_total = 0;
    int nMatchings1FTOF1_DifferentSectors_total = 0;
    int nMatchings1FTOF2_DifferentSectors_total = 0;

    int nMatchings2FTOF_total = 0;
    int nMatchings2FTOF1_total = 0;
    int nMatchings2FTOF2_total = 0;
    int nMatchings2FTOF_DifferentSectors_total = 0;
    int nMatchings2FTOF1_DifferentSectors_total = 0;
    int nMatchings2FTOF2_DifferentSectors_total = 0;

    int nMatchings3FTOF_total = 0;
    int nMatchings3FTOF1_total = 0;
    int nMatchings3FTOF2_total = 0;
    int nMatchings3FTOF_DifferentSectors_total = 0;
    int nMatchings3FTOF1_DifferentSectors_total = 0;
    int nMatchings3FTOF2_DifferentSectors_total = 0;

    /* Matching DC-ECAL */
    CrossMatcherToECalClusters ECMatcher;
    private double minDistance_ECAL = 50.;

    /* Matching DC-FTOF2 */
    CrossMatcherToReconFTOFHits FTOF2Matcher;
    private double minDistanceCSI_FTOF2 = 5.;
    private double minDistanceY_FTOF2 = 10.;
    private double minDistance_FTOF2 = 50.;

    /* Matching DC-FTOF1B */
    CrossMatcherToReconFTOFHits FTOF1BMatcher;
    private double minDistanceCSI_FTOF1B = 5.;
    private double minDistanceY_FTOF1B = 10.;
    private double minDistance_FTOF1B = 50.;

    /* Matching DC-FTOF1A */
    CrossMatcherToReconFTOFHits FTOF1AMatcher;
    private double minDistanceCSI_FTOF1A = 5.;
    private double minDistanceY_FTOF1A = 10.;
    private double minDistance_FTOF1A = 510.;

    /* Time history */
    private int TH_multiplicity;
    private double TH_minE_FTOF;
    private double TH_minE_PCAL;
    private double TH_deltaT;
    private double TH_deltaR;
    private int TH_DC;

    /* Momentum array */
    ArrayList<Double> Parray; // Index of lower-bound momenta (a part from
                              // latest, that is last bin
                              // higher-bound momentum)

    /* Variables */
    boolean doShowHistograms = false;
    int nevent = -1;

    int nTracksMatchedToGen = 1; /* To avoid divide-by-0 */
    int nTracksQPMatchedToGen = 1;
    int nTracksQMMatchedToGen = 1;

    int nEvents1Rec = 0;
    int nEvents2Rec = 0;
    int nEvents3Rec = 0;

    ArrayList<Integer> nTracksMatchedToGen_vsP = new ArrayList<Integer>();
    ArrayList<Integer> nTracksMatchedToGenQP_vsP = new ArrayList<Integer>();
    ArrayList<Integer> nTracksMatchedToGenQM_vsP = new ArrayList<Integer>();

    /* Histogram class */
    GuiClass guiClass;

    public static void main(String[] args) throws IOException {
        // TODO Auto-generated method stub
        AnalysisClass myClass = new AnalysisClass();
        myClass.setup(args);
        myClass.print();
        myClass.run();
    }

    /* Matching methods - working on the MatchedCrosses */
    private void matchToClusters(List<? extends MatchedCross> crosses, List<ECCluster> clusters[]) {
        double dMatching;
        for (MatchedCross cross : crosses) {
            dMatching = this.ECMatcher.matchCrossToHits(cross, clusters[cross.get_Sector() - 1]);
            cross.setMinDistanceEC(dMatching);
        }
    }

    private void matchToFTOFHits(int layer, List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {
        switch (layer) {
        case 1:
            this.matchToFTOF1AHits(crosses, hits);
            break;
        case 2:
            this.matchToFTOF1BHits(crosses, hits);
            break;
        case 3:
            this.matchToFTOF2Hits(crosses, hits);
            break;
        }
        return;
    }

    private void matchToFTOF1AHits(List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {
        double dist = -1;

        for (MatchedCross cross : crosses) {
            dist = this.FTOF1AMatcher.matchCrossToHits(cross, hits[cross.get_Sector() - 1]);
            if (this.FTOF1AMatcher.distanceIsSmallerThanMin(dist)) {
                cross.setIsMatchedToFTOF1A(true);
                cross.setMinDistanceFTOF1A(dist);
            }
        }
    }

    private void matchToFTOF1BHits(List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {
        double dist = -1;

        for (MatchedCross cross : crosses) {
            dist = this.FTOF1BMatcher.matchCrossToHits(cross, hits[cross.get_Sector() - 1]);
            if (this.FTOF1BMatcher.distanceIsSmallerThanMin(dist)) {
                cross.setIsMatchedToFTOF1B(true);
                cross.setMinDistanceFTOF1B(dist);
            }
        }
    }

    private void matchToFTOF2Hits(List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {

        double dist = -1;
        for (MatchedCross cross : crosses) {
            dist = this.FTOF2Matcher.matchCrossToHits(cross, hits[cross.get_Sector() - 1]);
            if (this.FTOF2Matcher.distanceIsSmallerThanMin(dist)) {
                cross.setIsMatchedToFTOF1B(true);
                cross.setMinDistanceFTOF1B(dist);
            }
        }
    }

    private void doTimeHistory() {
        nMatchingsTotalALL_thisEvent = 0;
        nMatchingsDiffSectorsALL_thisEvent = 0;

        nMatchingsTotalFTOF1_thisEvent = 0;
        nMatchingsDiffSectorsFTOF1_thisEvent = 0;

        nMatchingsTotalFTOF2_thisEvent = 0;
        nMatchingsDiffSectorsFTOF2_thisEvent = 0;

        /* Check FTOF1 / FTOF2 */
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            nFTOF1PCAL_thisEvent[ii] = 0;
            nFTOF2_thisEvent[ii] = 0;

            /* 1A-1B coinc. */
            List<HitsMatcherFTOF1PCal.MatchedFTOF1PCALHits> retList = this.hitsMatcherFTOF1[ii].MatchHits(hitsFTOF1A[ii], hitsFTOF1B[ii], clusters[ii]);
            for (HitsMatcherFTOF1PCal.MatchedFTOF1PCALHits hit : retList) {
                // System.out.println("TRIG --- -- "+nevent+" "+ii+" "+hit.getTime());
                if ((hit.getTime() > this.Tmin) && (hit.getTime() < this.Tmin + this.Tcoinc)) {
                    nFTOF1PCAL_thisEvent[ii]++;
                }
            }
            // System.out.println("RESULT TRIG SSUM: "+nevent+" "+ii+" "+matchFTOF1PCAL[ii]);
            /* 2 */
            for (RawTOFHit hit2 : hitsFTOF2[ii]) {        
                if ((hit2.get_TimeAVG() > this.Tmin) && (hit2.get_TimeAVG() < this.Tmin + this.Tcoinc) && (hit2.getEnergy()>this.TH_minE_FTOF)) {
                    if ((TH_DC == 0) || ((TH_DC == 1) && sectorHasDCsegments[ii]) || ((TH_DC >= 2) && hit2.isMatchedToR3CrossProjection())) nFTOF2_thisEvent[ii]++;
                }
            }

            nMatchingsTotalFTOF1_thisEvent += nFTOF1PCAL_thisEvent[ii];
            nMatchingsTotalFTOF2_thisEvent += nFTOF2_thisEvent[ii];
            nMatchingsTotalALL_thisEvent += (nFTOF2_thisEvent[ii] + nFTOF1PCAL_thisEvent[ii]);

            nMatchingsDiffSectorsFTOF1_thisEvent += (nFTOF1PCAL_thisEvent[ii] == 0) ? 0 : 1;
            nMatchingsDiffSectorsFTOF2_thisEvent += (nFTOF2_thisEvent[ii] == 0) ? 0 : 1;
            nMatchingsDiffSectorsALL_thisEvent += ((nFTOF2_thisEvent[ii] == 0) & (nFTOF1PCAL_thisEvent[ii] == 0)) ? 0 : 1;

        }
    }

    /* Direct analysis */
    /* Should be called only when there's one gen. particle only! */
    private void doDirectAnalysis() {

        if (genParticles.size() > 1) return;
        MatchedParticle particle = genParticles.get(0);
        double theta, phi, p;
        int imom;
        int q;
        int pid;

        int nCoincFTOF1 = 0;
        Integer tmpI;

        /* Get the variables */
        theta = Math.toDegrees(particle.theta());
        phi = Math.toDegrees(particle.phi());
        p = particle.p();
        q = particle.charge();
        pid = particle.pid();

        imom = this.getMomentumIndex(p);
        if (imom < 0) return;

        /* Fill the "generated histograms" */
        guiClass.getHistogram2D("h2_ThetaPhiAllGEN").fill(phi, theta);
        guiClass.getHistogram2D("h2_ThetaPhiAllGEN_" + imom).fill(phi, theta);

        if (q > 0) {
            guiClass.getHistogram2D("h2_ThetaPhiQPGEN").fill(phi, theta);
            guiClass.getHistogram2D("h2_ThetaPhiQPGEN_" + imom).fill(phi, theta);
        }
        if (q < 0) {
            guiClass.getHistogram2D("h2_ThetaPhiQMGEN").fill(phi, theta);
            guiClass.getHistogram2D("h2_ThetaPhiQMGEN_" + imom).fill(phi, theta);
        }
        /*
         * Now check if, among the reconstructed tracks, there's one matching to this gen. particle. If so, fill the "vs Distance" histograms, that are usefull
         * to decide the cross projection
         */
        for (TrackMatchedToGen track : tracksDC) {
            if (track.getIdGen() == 0) { /* We enter here is there's 1 gen. particle only, hence track should match this */
                nTracksMatchedToGen++;
                tmpI = nTracksMatchedToGen_vsP.get(imom);
                tmpI++;
                nTracksMatchedToGen_vsP.set(imom, tmpI);
                if (q > 0) {
                    nTracksQPMatchedToGen++;
                    tmpI = nTracksMatchedToGenQP_vsP.get(imom);
                    tmpI++;
                    nTracksMatchedToGenQP_vsP.set(imom, tmpI);
                }
                if (q < 0) {
                    nTracksQMMatchedToGen++;
                    tmpI = nTracksMatchedToGenQM_vsP.get(imom);
                    tmpI++;
                    nTracksMatchedToGenQM_vsP.set(imom, tmpI);
                }

                /* Fill the "vsDistance" histograms for the EC-matching */
                for (int ibin = 0; ibin < guiClass.getHistogram1D("h1_vsDistanceAllREC").getxAxis().getNBins(); ibin++) {
                    double thisD = guiClass.getHistogram1D("h1_vsDistanceAllREC").getxAxis().getBinCenter(ibin);
                    if (track.getMinDistanceEC() < thisD) {
                        if (q > 0) {
                            guiClass.getHistogram1D("h1_vsDistanceQPREC").incrementBinContent(ibin);
                            guiClass.getHistogram1D("h1_vsDistanceQPREC_" + imom).incrementBinContent(ibin);
                        }
                        if (q < 0) {
                            guiClass.getHistogram1D("h1_vsDistanceQMREC").incrementBinContent(ibin);
                            guiClass.getHistogram1D("h1_vsDistanceQMREC_" + imom).incrementBinContent(ibin);
                        }
                        guiClass.getHistogram1D("h1_vsDistanceAllREC").incrementBinContent(ibin);
                        guiClass.getHistogram1D("h1_vsDistanceAllREC_" + imom).incrementBinContent(ibin);
                    }
                }
                break;
            }
        }

        /*
         * Move forward only if this generated particle has been reconstructed
         */
        if (particle.isMatched()) {

           
            guiClass.getHistogram2D("h2_ThetaPhiAllREC").fill(phi, theta);
            guiClass.getHistogram2D("h2_ThetaPhiAllREC_" + imom).fill(phi, theta);

            if (q > 0) {
                guiClass.getHistogram2D("h2_ThetaPhiQPREC").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiQPREC_" + imom).fill(phi, theta);
            }
            if (q < 0) {
                guiClass.getHistogram2D("h2_ThetaPhiQMREC").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiQMREC_" + imom).fill(phi, theta);
            }

            if (nMatchingsTotalFTOF1_thisEvent >= 1) {

                guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER_" + imom).fill(phi, theta);

                if (q > 0) {
                    guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER_" + imom).fill(phi, theta);
                }
                if (q < 0) {
                    guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER_" + imom).fill(phi, theta);
                }
            }
            if (nMatchingsTotalALL_thisEvent >= 1) {
                guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER2").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER2_" + imom).fill(phi, theta);

                if (q > 0) {
                    guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER2").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER2_" + imom).fill(phi, theta);
                }
                if (q < 0) {
                    guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER2").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER2_" + imom).fill(phi, theta);
                }
            }

            else {
                 if ((p > 5) && (theta > 20) && (theta < 30) && (phi > -5) && (phi < 5)) {
                     System.out.println(nevent);
                 }
            }
        }/* end genParticleIsMatched (to a recon) */
    }

    /* Inverse analysis */
    private void doInverseCrossAnalysis() {

        ArrayList<Double> matchingDistances = new ArrayList<Double>();

        double dMin, thisD;
        int nMatchingsEC = 0;
        int nMatchingsFTOF2 = 0;

        /*
         * Want to check that the SAME cross is not in coincidence with FTOF2 and EC?
         */
        for (MatchedCross cross : crossesDC) {

            if (cross.isMatchedToFTOF2()) nMatchingsFTOF2++;

            dMin = cross.getMinDistanceEC();
            if (dMin < 0)
                continue; // No ECclusters at all have been found in this
                          // sector;
            else {
                matchingDistances.add(new Double(dMin));
            }
        }

        /*
         * Now, matchingDistances contains, for each cross in this event that has been matched, the corresponding distance. nMatchingsFTOF2 contains the number of crosses matching FTOF2
         */
        for (int ibin = 0; ibin < guiClass.getHistogram1D("h1_vsDistance1TRIGGER").getxAxis().getNBins(); ibin++) {
            thisD = guiClass.getHistogram1D("h1_vsDistance1TRIGGER").getxAxis().getBinCenter(ibin);
            nMatchingsEC = 0;
            for (Double matching : matchingDistances) {
                if (matching < thisD) nMatchingsEC++;
            }
            if (nMatchingsEC >= 1) guiClass.getHistogram1D("h1_vsDistance1TRIGGER").incrementBinContent(ibin);
            if (nMatchingsEC >= 2) guiClass.getHistogram1D("h1_vsDistance2TRIGGER").incrementBinContent(ibin);
            if (nMatchingsEC >= 3) guiClass.getHistogram1D("h1_vsDistance3TRIGGER").incrementBinContent(ibin);

            if ((nMatchingsEC + nMatchingsFTOF2) >= 1) guiClass.getHistogram1D("h1_vsDistance1TRIGGER2").incrementBinContent(ibin);
            if ((nMatchingsEC + nMatchingsFTOF2) >= 2) guiClass.getHistogram1D("h1_vsDistance2TRIGGER2").incrementBinContent(ibin);
            if ((nMatchingsEC + nMatchingsFTOF2) >= 3) guiClass.getHistogram1D("h1_vsDistance3TRIGGER2").incrementBinContent(ibin);
        }
    }

    private void doInverseTimeAnalysis() {

        /* FTOF1 only */
        if (nMatchingsTotalFTOF1_thisEvent >= 1) nMatchings1FTOF1_total++;
        if (nMatchingsTotalFTOF1_thisEvent >= 2) nMatchings2FTOF1_total++;
        if (nMatchingsTotalFTOF1_thisEvent >= 3) nMatchings3FTOF1_total++;

        if (nMatchingsDiffSectorsFTOF1_thisEvent >= 1) nMatchings1FTOF1_DifferentSectors_total++;
        if (nMatchingsDiffSectorsFTOF1_thisEvent >= 2) nMatchings2FTOF1_DifferentSectors_total++;
        if (nMatchingsDiffSectorsFTOF1_thisEvent >= 3) nMatchings3FTOF1_DifferentSectors_total++;

        /* FTOF2 only */
        if (nMatchingsTotalFTOF2_thisEvent >= 1) nMatchings1FTOF2_total++;
        if (nMatchingsTotalFTOF2_thisEvent >= 2) nMatchings2FTOF2_total++;
        if (nMatchingsTotalFTOF2_thisEvent >= 3) nMatchings3FTOF2_total++;

        if (nMatchingsDiffSectorsFTOF2_thisEvent >= 1) nMatchings1FTOF2_DifferentSectors_total++;
        if (nMatchingsDiffSectorsFTOF2_thisEvent >= 2) nMatchings2FTOF2_DifferentSectors_total++;
        if (nMatchingsDiffSectorsFTOF2_thisEvent >= 3) nMatchings3FTOF2_DifferentSectors_total++;

        /* All FTOF */
        if (nMatchingsTotalALL_thisEvent >= 1) nMatchings1FTOF_total++;
        if (nMatchingsTotalALL_thisEvent >= 2) nMatchings2FTOF_total++;
        if (nMatchingsTotalALL_thisEvent >= 3) nMatchings3FTOF_total++;

        if (nMatchingsDiffSectorsALL_thisEvent >= 1) nMatchings1FTOF_DifferentSectors_total++;
        if (nMatchingsDiffSectorsALL_thisEvent >= 2) nMatchings2FTOF_DifferentSectors_total++;
        if (nMatchingsDiffSectorsALL_thisEvent >= 3) nMatchings3FTOF_DifferentSectors_total++;

        
       
        
    }

    private void endTimeAnalysis() {
        double R1_FTOF1 = (nMatchings1FTOF1_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF1 = (nMatchings2FTOF1_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF1 = (nMatchings3FTOF1_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF1_DC = (nMatchings1FTOF1_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF1_DC = (nMatchings2FTOF1_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF1_DC = (nMatchings3FTOF1_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF2 = (nMatchings1FTOF2_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF2 = (nMatchings2FTOF2_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF2 = (nMatchings3FTOF2_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF2_DC = (nMatchings1FTOF2_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF2_DC = (nMatchings2FTOF2_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF2_DC = (nMatchings3FTOF2_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF = (nMatchings1FTOF_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF = (nMatchings2FTOF_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF = (nMatchings3FTOF_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF_DC = (nMatchings1FTOF_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF_DC = (nMatchings2FTOF_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF_DC = (nMatchings3FTOF_DifferentSectors_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        System.out.println("coinc. rate in CLAS12 FTOF1 all sectors (1 - 2  - 3 ch particles) kHz: " + R1_FTOF1 + " - " + R2_FTOF1 + " - " + R3_FTOF1);
        System.out.println("coinc. rate in CLAS12 FTOF1-Different sectors (1 - 2  - 3 ch particles) kHz: " + R1_FTOF1_DC + " - " + R2_FTOF1_DC + " - " + R3_FTOF1_DC);

        System.out.println("coinc. rate in CLAS12 FTOF2 all sectors (1 - 2  - 3 ch particles) kHz: " + R1_FTOF2 + " - " + R2_FTOF2 + " - " + R3_FTOF2);
        System.out.println("coinc. rate in CLAS12 FTOF2-Different sectors (1 - 2  - 3 ch particles) kHz: " + R1_FTOF2_DC + " - " + R2_FTOF2_DC + " - " + R3_FTOF2_DC);

        System.out.println("coinc. rate in CLAS12 FTOF all sectors (1 - 2  - 3 ch particles) kHz: " + R1_FTOF + " - " + R2_FTOF + " - " + R3_FTOF);
        System.out.println("coinc. rate in CLAS12 FTOF-Different sectors (1 - 2  - 3 ch particles) kHz: " + R1_FTOF_DC + " - " + R2_FTOF_DC + " - " + R3_FTOF_DC);

      
        
    }

    /* Here starts non-physics related methods */
    public void print() {
        System.out.println("Ana class");
        System.out.println("HIPO fname: " + this.inputFileName);
    }

    public void setup(String[] args) throws IOException {
        FileReader in = new FileReader(args[0]);
        BufferedReader br = new BufferedReader(in);
        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#")) /* This is a comment */
            continue;
            String[] splited = line.split(" ");
            if (splited[0].contains("fname")) {
                this.inputFileName = splited[1];
            } else if (splited[0].contains("nevents")) {
                this.nEvents = Integer.parseInt(splited[1]);
            } else if (splited[0].contains("showHistograms")) {
                this.doShowHistograms = true;
            } else if (splited[0].contains("idEC")) {
                this.idEC = Integer.parseInt(splited[1]);
            } else if (splited[0].contains("minClusterE_ECAL")) {
                this.minClusterE_ECAL = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistance_ECAL")) {
                this.minDistance_ECAL = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minE_FTOF2")) {
                this.minE_FTOF2 = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistanceCSI_FTOF2")) {
                this.minDistanceCSI_FTOF2 = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistanceY_FTOF2")) {
                this.minDistanceY_FTOF2 = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minE_FTOF1B")) {
                this.minE_FTOF1B = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistanceCSI_FTOF1B")) {
                this.minDistanceCSI_FTOF1B = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistanceY_FTOF1B")) {
                this.minDistanceY_FTOF1B = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minE_FTOF1A")) {
                this.minE_FTOF1A = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistanceCSI_FTOF1A")) {
                this.minDistanceCSI_FTOF1A = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistanceY_FTOF1A")) {
                this.minDistanceY_FTOF1A = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistance_FTOF1A")) {
                this.minDistance_FTOF1A = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistance_FTOF1B")) {
                this.minDistance_FTOF1B = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("minDistance_FTOF2")) {
                this.minDistance_FTOF2 = Double.parseDouble(splited[1]);
            }

            else if (splited[0].contains("Tmin")) {
                this.Tmin = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("Tcoinc")) {
                this.Tcoinc = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("TH_multiplicity")) {
                this.TH_multiplicity = Integer.parseInt(splited[1]);
            } else if (splited[0].contains("TH_minE_FTOF")) {
                this.TH_minE_FTOF = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("TH_minE_PCAL")) {
                this.TH_minE_PCAL = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("TH_deltaT")) {
                this.TH_deltaT = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("TH_deltaR")) {
                this.TH_deltaR = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("TH_DC")) {
                this.TH_DC = Integer.parseInt(splited[1]);
            }

        }

        br.close();
        this.setupMomentum();
        this.setupHistograms();
        this.setupGeo();
        this.setupDataReaderAndMatcher();

    }

    private void setupDataReaderAndMatcher() {
        dataReader = new DataReader(this);

        dataReader.setMinClusterE_ECAL(minClusterE_ECAL);
        dataReader.setMinE_FTOF2(minE_FTOF2);
        dataReader.setMinE_FTOF1B(minE_FTOF1B);
        dataReader.setMinE_FTOF1A(minE_FTOF1A);
        dataReader.setMinE_FTOF(this.TH_minE_FTOF);

        /* Also setup here what needed to read data */
        genParticles = new ArrayList<MatchedParticle>();
        recParticles = new ArrayList<MatchedParticle>();

        clusters = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            clusters[ii] = new ArrayList<ECCluster>();
        }

        tracksDC = new ArrayList<TrackMatchedToGen>();
        crossesDC = new ArrayList<MatchedCross>();
        segmentsDC = new ArrayList<SimpleDCSegment>();

        hitsFTOF2 = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF2[ii] = new ArrayList<ReconTOFHit>();
        }

        hitsFTOF1B = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF1B[ii] = new ArrayList<ReconTOFHit>();
        }
        hitsFTOF1A = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF1A[ii] = new ArrayList<ReconTOFHit>();
        }

        hitsMatcherFTOF1 = new HitsMatcherFTOF1PCal[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsMatcherFTOF1[ii] = new HitsMatcherFTOF1PCal(this, ii, this.TH_multiplicity, this.TH_deltaT, this.TH_deltaR, this.TH_minE_FTOF, this.TH_minE_PCAL, this.TH_DC);
            hitsMatcherFTOF1[ii].setRejectPCAL(true);
            hitsMatcherFTOF1[ii].setRejectFTOF1A(true);
        }

        sectorHasDCsegments = new boolean[6];

    }

    private void setupGeo() {

        this.setupGeoEC();
        this.setupGeoFTOF2();
        this.setupGeoFTOF1B();
        this.setupGeoFTOF1A();

    }

    private void setupGeoEC() {
        this.ECMatcher = new CrossMatcherToECalClusters(this, this.idEC);
        this.ECMatcher.setMinDistance(this.minDistance_ECAL);
        this.ECMatcher.setH1_matchedHitsE("h1_matchedPCAL_E");
        this.ECMatcher.setH1_closerDistance("h1_minPCAL_Distance");
        this.ECMatcher.setupGeo();
    }

    private void setupGeoFTOF2() {
        this.FTOF2Matcher = new CrossMatcherToReconFTOFHits(this, 3);
        // this.FTOF2Matcher.setMinDistanceCSI(this.minDistanceCSI_FTOF2);
        // this.FTOF2Matcher.setMinDistanceY(this.minDistanceY_FTOF2);
        this.FTOF2Matcher.setMinDistance(this.minDistance_FTOF2);
        this.FTOF2Matcher.setH1_matchedHitsE("h1_matchedTOF2_E");
        this.FTOF2Matcher.setH1_closerDistance("h1_minTOF2_Distance");
        this.FTOF2Matcher.setupGeo();
    }

    private void setupGeoFTOF1B() {
        this.FTOF1BMatcher = new CrossMatcherToReconFTOFHits(this, 2);
        // this.FTOF1BMatcher.setMinDistanceCSI(this.minDistanceCSI_FTOF1B);
        // this.FTOF1BMatcher.setMinDistanceY(this.minDistanceY_FTOF1B);
        this.FTOF1BMatcher.setH1_matchedHitsE("h1_matchedTOF1B_E");
        this.FTOF1BMatcher.setH1_closerDistance("h1_minTOF1B_Distance");
        this.FTOF1BMatcher.setMinDistance(this.minDistance_FTOF1B);
        this.FTOF1BMatcher.setupGeo();
    }

    private void setupGeoFTOF1A() {
        this.FTOF1AMatcher = new CrossMatcherToReconFTOFHits(this, 1);
        // this.FTOF1AMatcher.setMinDistanceCSI(this.minDistanceCSI_FTOF1A);
        // this.FTOF1AMatcher.setMinDistanceY(this.minDistanceY_FTOF1A);
        this.FTOF1AMatcher.setH1_matchedHitsE("h1_matchedTOF1A_E");
        this.FTOF1AMatcher.setH1_closerDistance("h1_minTOF1A_Distance");
        this.FTOF1AMatcher.setMinDistance(this.minDistance_FTOF1A);
        this.FTOF1AMatcher.setupGeo();
    }

    public H1F getHistogram1D(String name) {
        return guiClass.getHistogram1D(name);
    }

    public H2F getHistogram2D(String name) {
        return guiClass.getHistogram2D(name);
    }

    public void setupHistograms() {
        System.out.println("Going to setup histograms");
        guiClass = new GuiClass(this);
        guiClass.setParray(Parray);
        guiClass.setDeltaT(250.*1E-9);
        guiClass.setupHistograms();
        System.out.println("Setup histograms done");
    }

    public void setupMomentum() {
        /* Dirt and ugly - by hand */
        Parray = new ArrayList<Double>();
        Parray.add(new Double(0.3));
        Parray.add(new Double(1.));
        Parray.add(new Double(2.));
        Parray.add(new Double(3.));
        Parray.add(new Double(5.));
        Parray.add(new Double(11.));

        for (int ii = 0; ii < Parray.size() - 1; ii++) {
            nTracksMatchedToGen_vsP.add(new Integer(1)); /* to avoid divide-by-0 */
            nTracksMatchedToGenQP_vsP.add(new Integer(1));
            nTracksMatchedToGenQM_vsP.add(new Integer(1));
        }
    }

    public int getMomentumIndex(double P) {
        int index = -1;
        for (int ii = 0; ii < Parray.size() - 1; ii++) {
            if ((P >= Parray.get(ii)) && (P < Parray.get(ii + 1))) {
                index = ii;
                return index;
            }
        }
        return index;
    }

    public void run() {

        HipoDataSource reader = new HipoDataSource();
        reader.open(this.inputFileName);
        System.out.println("There are: " + reader.getSize() + " events in the file");
        H1F hi = new H1F("hi_ftof_energy_PA", "hi_ftof_energy_PA", 100, 0.0, 1.0);

        int nmultitrk = 0;

        if (nEvents == -1) {
            nEvents = reader.getSize();
            System.out.println("Reading ALL ");
        } else {
            System.out.println("Reading: " + nEvents);
        }

        while (reader.hasEvent() == true && nevent < nEvents) {

            nevent++;
            DataEvent event = reader.getNextEvent();
            if (nevent == 0) continue; // skip evt 0

            if (nevent % 10000 == 0) System.out.println("Analyzing " + nevent + " events" + " " + nmultitrk);

            /* Get generated (MC) particles */
            DataBank genParticlesDB = event.getBank("MC::Particle");
            DataBank recParticlesDB = event.getBank("REC::Particle");
            DataBank recParticlesTOFDB = event.getBank("REC::Scintillator");
            DataBank recParticlesCALDB = event.getBank("REC::Calorimeter");
            nGenParticles = dataReader.makeGeneratedParticles(genParticlesDB, genParticles);

            /* Get recon particles - those that ALSO have a FTOF hit (any) / a CALO hit (any) */
            nRecParticles = dataReader.makeReconstructedParticles(recParticlesDB, recParticlesTOFDB, recParticlesCALDB, recParticles);

            /* Match gen - rec */
            dataReader.matchReconstructedParticles(recParticles, genParticles);

            /* Make tracks and crosses and segments */
            DataBank tracksData = event.getBank("HitBasedTrkg::HBTracks");
            DataBank crossesData = event.getBank("HitBasedTrkg::HBCrosses");
            DataBank segmentsData = event.getBank("HitBasedTrkg::HBSegments");

            // DataBank tracksData = event.getBank("TimeBasedTrkg::TBTracks");
            // DataBank crossesData = event.getBank("TimeBasedTrkg::TBCrosses");

            tracksDC.clear();
            crossesDC.clear();
            segmentsDC.clear();

            dataReader.makeDCSegments(segmentsData, segmentsDC, sectorHasDCsegments);
            dataReader.makeDCCrossesAndTracks(tracksData, crossesData, tracksDC, crossesDC);

            /*
             * Check how many MC particles have been reconstructed and match generated.
             * This is used for the trigger efficiency: NtriggeredAndRecon / NRecon
             * -> I want to consider only good particles and check for them it they were triggered
             */
            int nReconstructedTracks = dataReader.matchReconstructedTracks(genParticles, tracksDC);

            if (tracksData.rows() >= 1) nEvents1Rec++;
            if (tracksData.rows() >= 2) nEvents2Rec++;
            if (tracksData.rows() >= 3) nEvents3Rec++;

            /* Get PCAL clusters */
            DataBank clustersEC = event.getBank("ECAL::clusters");
            int nClustersEC = dataReader.makeCaloClusters(idEC, clustersEC, clusters, TH_DC, sectorHasDCsegments);

            /*
             * Get FTOF hits - later will search for panel 2 (and maybe for panel 1A)
             */
            DataBank hitsRawFTOF = event.getBank("FTOF::rawhits");
            DataBank hitsReconFTOF = event.getBank("FTOF::hits");
            int nHitsFTOF2 = dataReader.makeFTOFHits(hitsRawFTOF, hitsReconFTOF, hitsFTOF2, hitsFTOF1B, hitsFTOF1A, TH_DC, sectorHasDCsegments);

            /* Perform matchings */         
            this.matchToClusters(crossesDC, clusters);
            this.matchToFTOFHits(3, crossesDC, hitsFTOF2); // FTOF2
            this.matchToFTOFHits(2, crossesDC, hitsFTOF1B); // FTOF1B
            this.matchToFTOFHits(1, crossesDC, hitsFTOF1A); // FTOF1A

            /*
             * Here do the "time-window analysis" that will determine trigger conditions.
             */ 
            this.doTimeHistory();

            /*
             * Do a "direct" analysis if nGen==1
             */
            if (nGenParticles == 1) {
                this.Tmin=0;
                this.Tcoinc=9999;
                this.doDirectAnalysis();
            }
            /*
             * Then, do an "inverse" analysis - basically what the trigger does: take all R3 crosses - doesn't matter if matched to a track or not, match them with ECAL.
             */

            this.doInverseCrossAnalysis();
            this.doInverseTimeAnalysis();

        } /* End loop on events */
        this.endTimeAnalysis();
        guiClass.processHistograms();
        if (this.doShowHistograms == true) {
            guiClass.showHistograms();
        }

    }

}
