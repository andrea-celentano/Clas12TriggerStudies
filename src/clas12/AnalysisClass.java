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
    List<Particle> genParticles;

    /* Tracks */
    List<TrackMatchedToGen> tracksDC;

    /* DC crosses */
    public static short crossIsAssociatedToTrack = 255;
    List<MatchedCross> crossesDC;

    /* DC segments */
    List<SimpleDCSegment> segmentsDC;

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

    int nMatchingsFTOF_thisEvent = 0; // 0 to 6
    int nMatchingsFTOF1_thisEvent = 0; // 0 to 6
    int nMatchingsFTOF2_thisEvent = 0; // 0 to 6

    int nMatchingsFTOF_DC_thisEvent = 0; // 0 to 6
    int nMatchingsFTOF1_DC_thisEvent = 0;
    int nMatchingsFTOF2_DC_thisEvent = 0; // 0 to 6

    int nMatchings1FTOF_total = 0;
    int nMatchings1FTOF1_total = 0;
    int nMatchings1FTOF2_total = 0;
    int nMatchings1FTOF_DC_total = 0;
    int nMatchings1FTOF1_DC_total = 0;
    int nMatchings1FTOF2_DC_total = 0;

    int nMatchings2FTOF_total = 0;
    int nMatchings2FTOF1_total = 0;
    int nMatchings2FTOF2_total = 0;
    int nMatchings2FTOF_DC_total = 0;
    int nMatchings2FTOF1_DC_total = 0;
    int nMatchings2FTOF2_DC_total = 0;

    int nMatchings3FTOF_total = 0;
    int nMatchings3FTOF1_total = 0;
    int nMatchings3FTOF2_total = 0;
    int nMatchings3FTOF_DC_total = 0;
    int nMatchings3FTOF1_DC_total = 0;
    int nMatchings3FTOF2_DC_total = 0;

    /* Matching DC-ECAL */
    CrossMatcherToECClusters ECMatcher;
    private double minDistance_ECAL = 50.;

    /* Matching DC-FTOF2 */
    CrossMatcherToRawFTOFHits FTOF2Matcher;
    private double minDistanceCSI_FTOF2 = 5.;
    private double minDistanceY_FTOF2 = 10.;

    /* Matching DC-FTOF1B */
    CrossMatcherToRawFTOFHits FTOF1BMatcher;
    private double minDistanceCSI_FTOF1B = 5.;
    private double minDistanceY_FTOF1B = 10.;

    /* Matching DC-FTOF1A */
    CrossMatcherToRawFTOFHits FTOF1AMatcher;
    private double minDistanceCSI_FTOF1A = 5.;
    private double minDistanceY_FTOF1A = 10.;

    /* Momentum array */
    ArrayList<Double> Parray; // Index of lower-bound momenta (a part from
                              // latest, that is last bin
                              // higher-bound momentum)

    /* Variables */
    boolean doShowHistograms = false;
    int nevent = -1;

    int nTracksWithR3Cross = 1; /* To avoid divide-by-0 */
    int nTracksQPWithR3Cross = 1;
    int nTracksQMWithR3Cross = 1;

    int nEvents1Rec = 0;
    int nEvents2Rec = 0;
    int nEvents3Rec = 0;

    ArrayList<Integer> nTracks_vsP = new ArrayList<Integer>();
    ArrayList<Integer> nTracksQP_vsP = new ArrayList<Integer>();
    ArrayList<Integer> nTracksQM_vsP = new ArrayList<Integer>();

    private int TH_multiplicity;
    private double TH_minE_FTOF;
    private double TH_minE_PCAL;
    private double TH_deltaT;
    private double TH_deltaR;

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
            dMatching = this.ECMatcher.matchCrossToClusters(cross, clusters[cross.get_Sector() - 1]);
            cross.setMinDistanceEC(dMatching);
        }
    }

    private void matchToFTOFHits(int layer, List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {
        switch (layer) {
        case 1:
        case 2:
            this.matchToFTOF1BHits(crosses, hits);
            return; // not yet implemented
        case 3:
            this.matchToFTOF2Hits(crosses, hits);
            return;
        }
    }

    private void matchToFTOF1BHits(List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {
        boolean isMatched = false;
        for (MatchedCross cross : crosses) {
            isMatched = this.FTOF1BMatcher.matchCrossToFTOFHits(cross, hits[cross.get_Sector() - 1]);
            cross.setIsMatchedToFTOF1B(isMatched);
        }
    }

    private void matchToFTOF2Hits(List<? extends MatchedCross> crosses, List<ReconTOFHit> hits[]) {
        boolean isMatched = false;
        for (MatchedCross cross : crosses) {
            isMatched = this.FTOF2Matcher.matchCrossToFTOFHits(cross, hits[cross.get_Sector() - 1]);
            cross.setIsMatchedToFTOF2(isMatched);
        }
    }

    private void doTimeHistory() {
        this.nMatchingsFTOF1_thisEvent = 0;
        this.nMatchingsFTOF2_thisEvent = 0;
        this.nMatchingsFTOF_thisEvent = 0;

        this.nMatchingsFTOF1_DC_thisEvent = 0;
        this.nMatchingsFTOF2_DC_thisEvent = 0;
        this.nMatchingsFTOF_DC_thisEvent = 0;

        int matchFTOF1PCAL[] = new int[AnalysisClass.nSectors_CLAS12];
        boolean hasFTOF2[] = new boolean[AnalysisClass.nSectors_CLAS12];
        boolean hasDCcrosses[] = new boolean[AnalysisClass.nSectors_CLAS12];
        boolean hasDCsegments[] = new boolean[AnalysisClass.nSectors_CLAS12];
        boolean hasDCsegmentsL1[] = new boolean[AnalysisClass.nSectors_CLAS12];
        boolean hasDCsegmentsL2[] = new boolean[AnalysisClass.nSectors_CLAS12];

        /*
         * Check if segments are in R3 - all segments, assuming DC are "slow" and time resolution (readout window)
         */
        for (SimpleDCSegment segment : segmentsDC) {
            if (segment.getLayer() == 3) {
                if (segment.getSuperlayer() == 5)
                    hasDCsegmentsL1[segment.getSector() - 1] = true;
                else if (segment.getSuperlayer() == 6) hasDCsegmentsL2[segment.getSector() - 1] = true;
            }
        }
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            if ((hasDCsegmentsL1[ii]) && (hasDCsegmentsL2[ii])) hasDCsegments[ii] = true;
        }

        for (Cross cross : crossesDC) {
            if (cross.get_Region() == 3) {
                hasDCcrosses[cross.get_Sector() - 1] = true;
            }
        }

        /* Check FTOF1 / FTOF2 */
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {

            /* 1A-1B coinc. */
            List<HitsMatcherFTOF1PCal.MatchedFTOF1PCALHits> retList = this.hitsMatcherFTOF1[ii].MatchHits(hitsFTOF1A[ii], hitsFTOF1B[ii], clusters[ii]);
            for (HitsMatcherFTOF1PCal.MatchedFTOF1PCALHits hit : retList) {
                // System.out.println("TRIG --- -- "+nevent+" "+ii+" "+hit.getTime());
                if ((hit.getTime() > this.Tmin) && (hit.getTime() < this.Tmin + this.Tcoinc)) {
                    matchFTOF1PCAL[ii]++;
                    break;
                }
            }
            // System.out.println("RESULT TRIG SSUM: "+nevent+" "+ii+" "+matchFTOF1PCAL[ii]);
            /* 2 */
            for (RawTOFHit hit2 : hitsFTOF2[ii]) {
                if ((hit2.get_TimeAVG() > this.Tmin) && (hit2.get_TimeAVG() < this.Tmin + this.Tcoinc)) {
                    hasFTOF2[ii] = true;
                    break;
                }
            }

            if (matchFTOF1PCAL[ii] > 0) this.nMatchingsFTOF1_thisEvent++;
            if ((matchFTOF1PCAL[ii] > 0) || hasFTOF2[ii]) this.nMatchingsFTOF_thisEvent++;
            if ((matchFTOF1PCAL[ii] > 0) && (hasDCcrosses[ii])) this.nMatchingsFTOF1_DC_thisEvent++;
            if (((matchFTOF1PCAL[ii] > 0) || hasFTOF2[ii]) && (hasDCcrosses[ii])) this.nMatchingsFTOF_DC_thisEvent++;

        }

    }

    /* Direct analysis */
    /* Should be called only when there's one gen. particle only! */
    private void doDirectAnalysis() {

        double theta, phi, p;
        int imom;
        int q;
        int pid;
        int idGen = -1;

        int nMatchingsCrossesEC = 0;
        int nMatchingsCrossesFTOF2 = 0;
        int nMatchingsCrossesFTOF1B = 0;

        int nCoincFTOF1 = 0;
        Integer tmpI;
        boolean genParticleIsRecon = false;

        for (Particle particle : genParticles) {

            idGen++;
            genParticleIsRecon = false;
            nMatchingsCrossesFTOF1B = 0;
            nMatchingsCrossesFTOF2 = 0;
            nMatchingsCrossesEC = 0;

            /* Get the variables */
            theta = Math.toDegrees(particle.theta());
            phi = Math.toDegrees(particle.phi());
            p = particle.p();
            q = particle.charge();
            pid = particle.pid();

            imom = this.getMomentumIndex(p);
            if (imom < 0) continue;

            /* Fill the "generated histograms" */
            guiClass.getHistogram2D("h2_ThetaPhiAllGEN").fill(phi, theta);
            guiClass.getHistogram2D("h2_ThetaPhiAllGEN_"+imom).fill(phi, theta);

            if (q > 0) {
                guiClass.getHistogram2D("h2_ThetaPhiQPGEN").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiQPGEN_"+imom).fill(phi, theta);
            }
            if (q < 0) {
                guiClass.getHistogram2D("h2_ThetaPhiQMGEN").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiQMGEN_"+imom).fill(phi, theta);
            }
            /*
             * Now check if, among the reconstructed tracks, there's one matching to this gen. particle. If so, fill the reconstructed histograms
             */
            for (TrackMatchedToGen track : tracksDC) {

                if (track.getIdGen() == idGen) {
                    genParticleIsRecon = true;

                    /* Also fill the "vsDistance" histograms for the EC-matching */
                    for (int ibin = 0; ibin <  guiClass.getHistogram1D("h1_vsDistanceAllREC").getxAxis().getNBins(); ibin++) {
                        double thisD = guiClass.getHistogram1D("h1_vsDistanceAllREC").getxAxis().getBinCenter(ibin);
                        if (track.getMinDistanceEC() < thisD) {
                            if (q > 0) {
                                guiClass.getHistogram1D("h1_vsDistanceQPREC").incrementBinContent(ibin);
                                guiClass.getHistogram1D("h1_vsDistanceQPREC_"+imom).incrementBinContent(ibin);
                            }
                            if (q < 0) {
                                guiClass.getHistogram1D("h1_vsDistanceQMREC").incrementBinContent(ibin);
                                guiClass.getHistogram1D("h1_vsDistanceQMREC_"+imom).incrementBinContent(ibin);
                            }
                            guiClass.getHistogram1D("h1_vsDistanceAllREC").incrementBinContent(ibin);
                            guiClass.getHistogram1D("h1_vsDistanceAllREC_"+imom).incrementBinContent(ibin);
                        }
                    }
                    break;
                }
            }

            /*
             * Move forward only if this generated particle has been reconstructed
             */
            if (genParticleIsRecon) {
                nTracksWithR3Cross++;
                tmpI = nTracks_vsP.get(imom);
                tmpI++;
                nTracks_vsP.set(imom, tmpI);
                guiClass.getHistogram2D("h2_ThetaPhiAllREC").fill(phi, theta);
                guiClass.getHistogram2D("h2_ThetaPhiAllREC_"+imom).fill(phi, theta);

                if (q > 0) {
                    nTracksQPWithR3Cross++;
                    tmpI = nTracksQP_vsP.get(imom);
                    tmpI++;
                    nTracksQP_vsP.set(imom, tmpI);

                    guiClass.getHistogram2D("h2_ThetaPhiQPREC").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiQPREC_"+imom).fill(phi, theta);
                }
                if (q < 0) {
                    nTracksQMWithR3Cross++;
                    tmpI = nTracksQM_vsP.get(imom);
                    tmpI++;
                    nTracksQM_vsP.set(imom, tmpI);

                    guiClass.getHistogram2D("h2_ThetaPhiQMREC").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiQMREC_"+imom).fill(phi, theta);
                }

                /* now count the number of cross-matchings */
                for (MatchedCross cross : crossesDC) {
                    if (this.ECMatcher.distanceIsSmallerThanMin(cross.getMinDistanceEC())) nMatchingsCrossesEC++;
                    if (cross.isMatchedToFTOF2()) nMatchingsCrossesFTOF2++;
                    if (cross.isMatchedToFTOF1B()) nMatchingsCrossesFTOF1B++;
                }

                if (nMatchingsCrossesFTOF1B >= genParticles.size()) {

                    guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER_"+imom).fill(phi, theta);

                    if (q > 0) {
                        guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER").fill(phi, theta);
                        guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER_"+imom).fill(phi, theta);
                    }
                    if (q < 0) {
                        guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER").fill(phi, theta);
                        guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER_"+imom).fill(phi, theta);
                    }
                }
                // if ((nMatchingsCrossesFTOF1B + nMatchingsCrossesFTOF2) >=
                // genParticles.size()) {
                if (nMatchingsFTOF_DC_thisEvent >= genParticles.size()) {
                    guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER2").fill(phi, theta);
                    guiClass.getHistogram2D("h2_ThetaPhiAllTRIGGER2_"+imom).fill(phi, theta);

                    if (q > 0) {
                        guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER2").fill(phi, theta);
                        guiClass.getHistogram2D("h2_ThetaPhiQPTRIGGER2_"+imom).fill(phi, theta);
                    }
                    if (q < 0) {
                        guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER2").fill(phi, theta);
                        guiClass.getHistogram2D("h2_ThetaPhiQMTRIGGER2_"+imom).fill(phi, theta);
                    }
                }

                else {
                    if ((p > 5) && (theta > 20) && (theta < 30) && (phi > 20) && (phi < 30)) {
                        System.out.println(nevent);
                    }
                }

            }/* end genParticleIsRecon */

        } /* end loop on genParticles */
    }

    /* Inverse analysis */
    private void doInverseCrossAnalysis() {

        this.matchToClusters(crossesDC, clusters);
        this.matchToFTOFHits(3, crossesDC, hitsFTOF2);
        this.matchToFTOFHits(2, crossesDC, hitsFTOF1B);

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
        for (int ibin = 0; ibin <   guiClass.getHistogram1D("h1_vsDistance1TRIGGER").getxAxis().getNBins(); ibin++) {
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
        if (nMatchingsFTOF1_thisEvent >= 1) nMatchings1FTOF1_total++;
        if (nMatchingsFTOF1_thisEvent >= 2) nMatchings2FTOF1_total++;
        if (nMatchingsFTOF1_thisEvent >= 3) nMatchings3FTOF1_total++;

        if (nMatchingsFTOF1_DC_thisEvent >= 1) nMatchings1FTOF1_DC_total++;
        if (nMatchingsFTOF1_DC_thisEvent >= 2) nMatchings2FTOF1_DC_total++;
        if (nMatchingsFTOF1_DC_thisEvent >= 3) nMatchings3FTOF1_DC_total++;

        /* All FTOF */
        if (nMatchingsFTOF_thisEvent >= 1) nMatchings1FTOF_total++;
        if (nMatchingsFTOF_thisEvent >= 2) nMatchings2FTOF_total++;
        if (nMatchingsFTOF_thisEvent >= 3) nMatchings3FTOF_total++;

        if (nMatchingsFTOF_DC_thisEvent >= 1) nMatchings1FTOF_DC_total++;
        if (nMatchingsFTOF_DC_thisEvent >= 2) nMatchings2FTOF_DC_total++;
        if (nMatchingsFTOF_DC_thisEvent >= 3) nMatchings3FTOF_DC_total++;

    }

    private void endTimeAnalysis() {
        double R1_FTOF1 = (nMatchings1FTOF1_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF1 = (nMatchings2FTOF1_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF1 = (nMatchings3FTOF1_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF1_DC = (nMatchings1FTOF1_DC_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF1_DC = (nMatchings2FTOF1_DC_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF1_DC = (nMatchings3FTOF1_DC_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF = (nMatchings1FTOF_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF = (nMatchings2FTOF_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF = (nMatchings3FTOF_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        double R1_FTOF_DC = (nMatchings1FTOF_DC_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R2_FTOF_DC = (nMatchings2FTOF_DC_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;
        double R3_FTOF_DC = (nMatchings3FTOF_DC_total / (this.nevent * this.Tcoinc * 1E-9)) * 1E-3;

        System.out.println("coinc. rate in CLAS12 FTOF1 only (1 - 2  - 3 ch particles) kHz: " + R1_FTOF1 + " - " + R2_FTOF1 + " - " + R3_FTOF1);
        System.out.println("coinc. rate in CLAS12 FTOF1-DC (1 - 2  - 3 ch particles) kHz: " + R1_FTOF1_DC + " - " + R2_FTOF1_DC + " - " + R3_FTOF1_DC);

        System.out.println("coinc. rate in CLAS12 FTOF only (1 - 2  - 3 ch particles) kHz: " + R1_FTOF + " - " + R2_FTOF + " - " + R3_FTOF);
        System.out.println("coinc. rate in CLAS12 FTOF-DC (1 - 2  - 3 ch particles) kHz: " + R1_FTOF_DC + " - " + R2_FTOF_DC + " - " + R3_FTOF_DC);

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
            } else if (splited[0].contains("Tmin")) {
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

        /* Also setup here what needed to read data */
        genParticles = new ArrayList<Particle>();

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
            hitsMatcherFTOF1[ii] = new HitsMatcherFTOF1PCal(this, ii, this.TH_multiplicity, this.TH_deltaT, this.TH_deltaR, this.TH_minE_FTOF, this.TH_minE_PCAL);
            hitsMatcherFTOF1[ii].setRejectPCAL(true);
        }

    }

    private void setupGeo() {

        this.setupGeoEC();
        this.setupGeoFTOF2();
        this.setupGeoFTOF1B();
        this.setupGeoFTOF1A();

    }

    private void setupGeoEC() {

        this.ECMatcher = new CrossMatcherToECClusters(this, this.idEC);
        this.ECMatcher.setupGeo();
        this.ECMatcher.setMinDistance(this.minDistance_ECAL);

    }

    private void setupGeoFTOF2() {

        this.FTOF2Matcher = new CrossMatcherToRawFTOFHits(this, 3);
        this.FTOF2Matcher.setupGeo();
        this.FTOF2Matcher.setMinDistanceCSI(this.minDistanceCSI_FTOF2);
        this.FTOF2Matcher.setMinDistanceY(this.minDistanceY_FTOF2);

    }

    private void setupGeoFTOF1B() {

        this.FTOF1BMatcher = new CrossMatcherToRawFTOFHits(this, 2);
        this.FTOF1BMatcher.setupGeo();
        this.FTOF1BMatcher.setMinDistanceCSI(this.minDistanceCSI_FTOF1B);
        this.FTOF1BMatcher.setMinDistanceY(this.minDistanceY_FTOF1B);

    }

    private void setupGeoFTOF1A() {

        this.FTOF1AMatcher = new CrossMatcherToRawFTOFHits(this, 1);
        this.FTOF1AMatcher.setupGeo();
        this.FTOF1AMatcher.setMinDistanceCSI(this.minDistanceCSI_FTOF1A);
        this.FTOF1AMatcher.setMinDistanceY(this.minDistanceY_FTOF1A);

    }

    
    public H1F getHistogram1D(String name) {     
        return guiClass.getHistogram1D(name);
    }

    public H2F getHistogram2D(String name) {
        return guiClass.getHistogram2D(name);
    }

    public void setupHistograms() {
        guiClass = new GuiClass(this);
        guiClass.setParray(Parray);
        guiClass.setupHistograms();
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
            nTracks_vsP.add(new Integer(1)); /* to avoid divide-by-0 */
            nTracksQP_vsP.add(new Integer(1));
            nTracksQM_vsP.add(new Integer(1));
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
            nGenParticles = dataReader.makeGeneratedParticles(genParticlesDB, genParticles);

            /* Make tracks and crosses and segments */
            DataBank tracksData = event.getBank("HitBasedTrkg::HBTracks");
            DataBank crossesData = event.getBank("HitBasedTrkg::HBCrosses");
            DataBank segmentsData = event.getBank("HitBasedTrkg::HBSegments");

            // DataBank tracksData = event.getBank("TimeBasedTrkg::TBTracks");
            // DataBank crossesData = event.getBank("TimeBasedTrkg::TBCrosses");

            tracksDC.clear();
            crossesDC.clear();
            segmentsDC.clear();

            dataReader.makeDCSegments(segmentsData, segmentsDC);
            dataReader.makeDCCrossesAndTracks(tracksData, crossesData, tracksDC, crossesDC);

            /*
             * Check how many MC particles have been reconstructed and match MC - Track
             */
            int nReconstructedParticles = dataReader.matchReconstructedTracks(genParticles, tracksDC);

            if (tracksData.rows() >= 1) nEvents1Rec++;
            if (tracksData.rows() >= 2) nEvents2Rec++;
            if (tracksData.rows() >= 3) nEvents3Rec++;

            /* Get calorimeter clusters - later will search for PCAL */
            DataBank clustersEC = event.getBank("ECAL::clusters");

            int nClustersEC = dataReader.makeCaloClusters(idEC, clustersEC, clusters);

            /*
             * Get FTOF hits - later will search for panel 2 (and maybe for panel 1A)
             */
            DataBank hitsRawFTOF = event.getBank("FTOF::rawhits");
            DataBank hitsReconFTOF = event.getBank("FTOF::hits");
            int nHitsFTOF2 = dataReader.makeFTOFHits(hitsRawFTOF, hitsReconFTOF, hitsFTOF2, hitsFTOF1B, hitsFTOF1A);

            /* Perform matchings */
            /* Crosses */
            this.matchToClusters(crossesDC, clusters);
            this.matchToFTOFHits(3, crossesDC, hitsFTOF2);
            this.matchToFTOFHits(2, crossesDC, hitsFTOF1B);

            /* Tracks */
            this.matchToClusters(tracksDC, clusters);
            this.matchToFTOFHits(3, tracksDC, hitsFTOF2);
            this.matchToFTOFHits(2, tracksDC, hitsFTOF1B);
            /*
             * Here do the "time-window analysis".
             */
            this.doTimeHistory();

            /*
             * First, do a "direct" analysis
             */
            if (nGenParticles == 1) {
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
