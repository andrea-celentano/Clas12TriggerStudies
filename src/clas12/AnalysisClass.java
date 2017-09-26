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

    private DataReaderAndMatcher dataReader;

    /* Here goes common definitions */
    public static double phiAngle_CLAS12 = 60.;
    public static double thetaAngleDC_CLAS12 = 25.;
    public static int nSectors_CLAS12 = 6;

    /* Generated particles */
    int nGenParticles;
    List<Particle> genParticles;

    /* Tracks */
    List<TrackMatchedToGen> tracks;

    /* Crosses */
    public static short crossIsAssociatedToTrack = 255;
    List<MatchedCross> crosses;

    /* Clusters */
    private int idEC = 0; // 0: PCAL, 1: EC-IN, 2: EC-OUT
    private double minClusterE_ECAL = 0.;
    List<ECCluster>[] clusters;

    /* TOF */
    private double minE_FTOF2 = 0.1;
    private double minE_FTOF1B = 0.1;
    private double minE_FTOF1A = 0.1;
    List<SimpleTOFHit>[] hitsFTOF2;
    List<SimpleTOFHit>[] hitsFTOF1B;
    List<SimpleTOFHit>[] hitsFTOF1A;

    /* TOF "time history matching */
    double Tmin=100;
    double Tcoinc=100;
    HitsMatcherFTOF1 hitsMatcherFTOF1[];
    boolean          matchFTOF1result[];
    /* Matching DC-ECAL */
    CrossMatcherToECClusters ECMatcher;
    private double minDistance_ECAL = 50.;

    /* Matching DC-FTOF2 */
    CrossMatcherToFTOFHits FTOF2Matcher;
    private double minDistanceCSI_FTOF2 = 5.;
    private double minDistanceY_FTOF2 = 10.;

    /* Matching DC-FTOF1B */
    CrossMatcherToFTOFHits FTOF1BMatcher;
    private double minDistanceCSI_FTOF1B = 5.;
    private double minDistanceY_FTOF1B = 10.;

    /* Matching DC-FTOF1A */
    CrossMatcherToFTOFHits FTOF1AMatcher;
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

    /* Here goes histograms */
    ArrayList<H1F> h1_minCrossClusterDistance;
    ArrayList<H1F> h1_allClustersE;
    ArrayList<H1F> h1_closerClustersE;

    H2F h2_ThetaPhiAllGEN;
    H2F h2_ThetaPhiAllREC;
    H2F h2_ThetaPhiAllTRIGGER;
    H2F h2_ThetaPhiAllTRIGGER2;
    H2F h2_ThetaPhiAllEFF;
    H2F h2_ThetaPhiAllEFF2;

    H2F h2_ThetaPhiQPGEN;
    H2F h2_ThetaPhiQPREC;
    H2F h2_ThetaPhiQPTRIGGER;
    H2F h2_ThetaPhiQPTRIGGER2;
    H2F h2_ThetaPhiQPEFF;
    H2F h2_ThetaPhiQPEFF2;

    H2F h2_ThetaPhiQMGEN;
    H2F h2_ThetaPhiQMREC;
    H2F h2_ThetaPhiQMTRIGGER;
    H2F h2_ThetaPhiQMTRIGGER2;
    H2F h2_ThetaPhiQMEFF;
    H2F h2_ThetaPhiQMEFF2;

    ArrayList<H2F> h2_ThetaPhiAllGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiAllREC_vsP;
    ArrayList<H2F> h2_ThetaPhiAllTRIGGER_vsP;
    ArrayList<H2F> h2_ThetaPhiAllTRIGGER2_vsP;
    ArrayList<H2F> h2_ThetaPhiAllEFF_vsP;
    ArrayList<H2F> h2_ThetaPhiAllEFF2_vsP;

    ArrayList<H2F> h2_ThetaPhiQPGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiQPREC_vsP;
    ArrayList<H2F> h2_ThetaPhiQPTRIGGER_vsP;
    ArrayList<H2F> h2_ThetaPhiQPTRIGGER2_vsP;
    ArrayList<H2F> h2_ThetaPhiQPEFF_vsP;
    ArrayList<H2F> h2_ThetaPhiQPEFF2_vsP;

    ArrayList<H2F> h2_ThetaPhiQMGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiQMREC_vsP;
    ArrayList<H2F> h2_ThetaPhiQMTRIGGER_vsP;
    ArrayList<H2F> h2_ThetaPhiQMTRIGGER2_vsP;
    ArrayList<H2F> h2_ThetaPhiQMEFF_vsP;
    ArrayList<H2F> h2_ThetaPhiQMEFF2_vsP;

    H1F h1_vsDistanceAllREC;
    H1F h1_vsDistanceAllTRIGGER;
    H1F h1_vsDistanceAllTRIGGER2;
    H1F h1_vsDistanceAllEFF;
    H1F h1_vsDistanceAllEFF2;

    H1F h1_vsDistanceQPREC;
    H1F h1_vsDistanceQPTRIGGER;
    H1F h1_vsDistanceQPTRIGGER2;
    H1F h1_vsDistanceQPEFF;
    H1F h1_vsDistanceQPEFF2;

    H1F h1_vsDistanceQMREC;
    H1F h1_vsDistanceQMTRIGGER;
    H1F h1_vsDistanceQMTRIGGER2;
    H1F h1_vsDistanceQMEFF;
    H1F h1_vsDistanceQMEFF2;

    ArrayList<H1F> h1_vsDistanceAllREC_vsP;
    ArrayList<H1F> h1_vsDistanceAllTRIGGER_vsP;
    ArrayList<H1F> h1_vsDistanceAllTRIGGER2_vsP;
    ArrayList<H1F> h1_vsDistanceAllEFF_vsP;
    ArrayList<H1F> h1_vsDistanceAllEFF2_vsP;

    ArrayList<H1F> h1_vsDistanceQPREC_vsP;
    ArrayList<H1F> h1_vsDistanceQPTRIGGER_vsP;
    ArrayList<H1F> h1_vsDistanceQPTRIGGER2_vsP;
    ArrayList<H1F> h1_vsDistanceQPEFF_vsP;
    ArrayList<H1F> h1_vsDistanceQPEFF2_vsP;

    ArrayList<H1F> h1_vsDistanceQMREC_vsP;
    ArrayList<H1F> h1_vsDistanceQMTRIGGER_vsP;
    ArrayList<H1F> h1_vsDistanceQMTRIGGER2_vsP;
    ArrayList<H1F> h1_vsDistanceQMEFF_vsP;
    ArrayList<H1F> h1_vsDistanceQMEFF2_vsP;

    H1F h1_vsDistance1TRIGGER;
    H1F h1_vsDistance1TRIGGER2;
    H1F h1_vsDistance1EFF;
    H1F h1_vsDistance1EFF2;

    H1F h1_vsDistance2TRIGGER;
    H1F h1_vsDistance2TRIGGER2;
    H1F h1_vsDistance2EFF;
    H1F h1_vsDistance2EFF2;

    H1F h1_vsDistance3TRIGGER;
    H1F h1_vsDistance3TRIGGER2;
    H1F h1_vsDistance3EFF;
    H1F h1_vsDistance3EFF2;

    H2F h2_FTOF2EnergyAll_LR;
    H2F h2_FTOF2EnergyMatched_LR;

    H2F h2_FTOF1BEnergyAll_LR;
    H2F h2_FTOF1BEnergyMatched_LR;

    H2F h2_FTOF1AEnergyAll_LR;
    H2F h2_FTOF1AEnergyMatched_LR;

    H2F h2tmp;

    ArrayList<H1F> allH1F;
    ArrayList<H2F> allH2F;

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

    private void matchToFTOFHits(int layer, List<? extends MatchedCross> crosses, List<SimpleTOFHit> hits[]) {
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

    private void matchToFTOF1BHits(List<? extends MatchedCross> crosses, List<SimpleTOFHit> hits[]) {
        boolean isMatched = false;
        for (MatchedCross cross : crosses) {
            isMatched = this.FTOF1BMatcher.matchCrossToFTOFHits(cross, hits[cross.get_Sector() - 1]);
            cross.setIsMatchedToFTOF1B(isMatched);
        }
    }

    private void matchToFTOF2Hits(List<? extends MatchedCross> crosses, List<SimpleTOFHit> hits[]) {
        boolean isMatched = false;
        for (MatchedCross cross : crosses) {
            isMatched = this.FTOF2Matcher.matchCrossToFTOFHits(cross, hits[cross.get_Sector() - 1]);
            cross.setIsMatchedToFTOF2(isMatched);
        }
    }

    /* Direct analysis */
    /* Should be called only when there's one gen. particle only! */

    /*
     * Fill following histograms (theta/phi - p bins): 1. Generated events 2.
     * Events with reconstructed track 3. Events that have (nMatchingCrossesEC
     * >= nGeneratedEvents) 4. Events that have
     * (nMatchingCrossesEC+nMatchingCrossesFTOF >=nGeneratedEvents)
     */
    private void doDirectAnalysis() {

        double theta, phi, p;
        int imom;
        int q;
        int pid;
        int idGen = 0;

        int nMatchingsCrossesEC = 0;
        int nMatchingsCrossesFTOF2 = 0;
        int nMatchingsCrossesFTOF1B = 0;

        Integer tmpI;

        for (Particle particle : genParticles) {

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
            h2_ThetaPhiAllGEN.fill(phi, theta);
            h2_ThetaPhiAllGEN_vsP.get(imom).fill(phi, theta);

            if (q > 0) {
                h2_ThetaPhiQPGEN.fill(phi, theta);
                h2_ThetaPhiQPGEN_vsP.get(imom).fill(phi, theta);
            }
            if (q < 0) {
                h2_ThetaPhiQMGEN.fill(phi, theta);
                h2_ThetaPhiQMGEN_vsP.get(imom).fill(phi, theta);
            }
            /*
             * Now check if, among the reconstructed tracks, there's one
             * matching to this gen. particle. If so, fill the reconstructed
             * histograms
             */
            for (TrackMatchedToGen track : tracks) {

                // if (track.getIdGen() == idGen) {
                if (true) {
                    nTracksWithR3Cross++;
                    tmpI = nTracks_vsP.get(imom);
                    tmpI++;
                    nTracks_vsP.set(imom, tmpI);
                    h2_ThetaPhiAllREC.fill(phi, theta);
                    h2_ThetaPhiAllREC_vsP.get(imom).fill(phi, theta);

                    if (q > 0) {
                        nTracksQPWithR3Cross++;
                        tmpI = nTracksQP_vsP.get(imom);
                        tmpI++;
                        nTracksQP_vsP.set(imom, tmpI);

                        h2_ThetaPhiQPREC.fill(phi, theta);
                        h2_ThetaPhiQPREC_vsP.get(imom).fill(phi, theta);
                    }
                    if (q < 0) {
                        nTracksQMWithR3Cross++;
                        tmpI = nTracksQM_vsP.get(imom);
                        tmpI++;
                        nTracksQM_vsP.set(imom, tmpI);

                        h2_ThetaPhiQMREC.fill(phi, theta);
                        h2_ThetaPhiQMREC_vsP.get(imom).fill(phi, theta);
                    }

                    /* Also fill the "vsDistance" histograms for the EC-matching */
                    for (int ibin = 0; ibin < h1_vsDistanceAllREC.getxAxis().getNBins(); ibin++) {
                        double thisD = h1_vsDistanceAllREC.getxAxis().getBinCenter(ibin);
                        if (track.getMinDistanceEC() < thisD) {
                            if (q > 0) {
                                h1_vsDistanceQPREC.incrementBinContent(ibin);
                                h1_vsDistanceQPREC_vsP.get(imom).incrementBinContent(ibin);
                            }
                            if (q < 0) {
                                h1_vsDistanceQMREC.incrementBinContent(ibin);
                                h1_vsDistanceQMREC_vsP.get(imom).incrementBinContent(ibin);
                            }
                            h1_vsDistanceAllREC.incrementBinContent(ibin);
                            h1_vsDistanceAllREC_vsP.get(imom).incrementBinContent(ibin);
                        }
                    }
                }
            }

            /* now count the number of cross-matchings */
            for (MatchedCross cross : crosses) {
                if (this.ECMatcher.distanceIsSmallerThanMin(cross.getMinDistanceEC())) nMatchingsCrossesEC++;
                if (cross.isMatchedToFTOF2()) nMatchingsCrossesFTOF2++;
                if (cross.isMatchedToFTOF1B()) nMatchingsCrossesFTOF1B++;

            }

            if (nMatchingsCrossesFTOF1B >= genParticles.size()) {

                h2_ThetaPhiAllTRIGGER.fill(phi, theta);
                h2_ThetaPhiAllTRIGGER_vsP.get(imom).fill(phi, theta);

                if (q > 0) {
                    h2_ThetaPhiQPTRIGGER.fill(phi, theta);
                    h2_ThetaPhiQPTRIGGER_vsP.get(imom).fill(phi, theta);
                }
                if (q < 0) {
                    h2_ThetaPhiQMTRIGGER.fill(phi, theta);
                    h2_ThetaPhiQMTRIGGER_vsP.get(imom).fill(phi, theta);
                }
            }
            if ((nMatchingsCrossesFTOF1B + nMatchingsCrossesFTOF2) >= genParticles.size()) {
                h2_ThetaPhiAllTRIGGER2.fill(phi, theta);
                h2_ThetaPhiAllTRIGGER2_vsP.get(imom).fill(phi, theta);

                if (q > 0) {
                    h2_ThetaPhiQPTRIGGER2.fill(phi, theta);
                    h2_ThetaPhiQPTRIGGER2_vsP.get(imom).fill(phi, theta);
                }
                if (q < 0) {
                    h2_ThetaPhiQMTRIGGER2.fill(phi, theta);
                    h2_ThetaPhiQMTRIGGER2_vsP.get(imom).fill(phi, theta);
                }

            }

            idGen++;
        }/* end loop on genParticles */

    }

    /* Inverse analysis */
    private void doInverseCrossAnalysis() {

        this.matchToClusters(crosses, clusters);
        this.matchToFTOFHits(3, crosses, hitsFTOF2);
        this.matchToFTOFHits(2, crosses, hitsFTOF1B);

        ArrayList<Double> matchingDistances = new ArrayList<Double>();

        double dMin, thisD;
        int nMatchingsEC = 0;
        int nMatchingsFTOF2 = 0;

        /*
         * Want to check that the SAME cross is not in coincidence with FTOF2
         * and EC?
         */
        for (MatchedCross cross : crosses) {

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
         * Now, matchingDistances contains, for each cross in this event that
         * has been matched, the corresponding distance. nMatchingsFTOF2
         * contains the number of crosses matching FTOF2
         */

        for (int ibin = 0; ibin < h1_vsDistance1TRIGGER.getxAxis().getNBins(); ibin++) {
            thisD = h1_vsDistance1TRIGGER.getxAxis().getBinCenter(ibin);
            nMatchingsEC = 0;
            for (Double matching : matchingDistances) {
                if (matching < thisD) nMatchingsEC++;
            }
            if (nMatchingsEC >= 1) h1_vsDistance1TRIGGER.incrementBinContent(ibin);
            if (nMatchingsEC >= 2) h1_vsDistance2TRIGGER.incrementBinContent(ibin);
            if (nMatchingsEC >= 3) h1_vsDistance3TRIGGER.incrementBinContent(ibin);

            if ((nMatchingsEC + nMatchingsFTOF2) >= 1) h1_vsDistance1TRIGGER2.incrementBinContent(ibin);
            if ((nMatchingsEC + nMatchingsFTOF2) >= 2) h1_vsDistance2TRIGGER2.incrementBinContent(ibin);
            if ((nMatchingsEC + nMatchingsFTOF2) >= 3) h1_vsDistance3TRIGGER2.incrementBinContent(ibin);

            if ((nMatchingsEC + nMatchingsFTOF2) >= 3) {
                System.out.println(nevent);
            }
        }
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
            }else if (splited[0].contains("Tmin")) {
                this.Tmin = Double.parseDouble(splited[1]);
            } else if (splited[0].contains("Tcoinc")) {
                this.Tcoinc = Double.parseDouble(splited[1]);
            } 

        }

        br.close();
        this.setupMomentum();
        this.setupHistograms();
        this.setupGeo();
        this.setupDataReaderAndMatcher();

    }

    private void setupDataReaderAndMatcher() {
        dataReader = new DataReaderAndMatcher(this);
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

        tracks = new ArrayList<TrackMatchedToGen>();
        crosses = new ArrayList<MatchedCross>();

        hitsFTOF2 = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF2[ii] = new ArrayList<SimpleTOFHit>();
        }

        hitsFTOF1B = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF1B[ii] = new ArrayList<SimpleTOFHit>();
        }
        hitsFTOF1A = new ArrayList[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF1A[ii] = new ArrayList<SimpleTOFHit>();
        }
        hitsMatcherFTOF1=new HitsMatcherFTOF1[AnalysisClass.nSectors_CLAS12];
        matchFTOF1result=new boolean[AnalysisClass.nSectors_CLAS12];
        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsMatcherFTOF1[ii]=new HitsMatcherFTOF1(this, ii, Tcoinc, Tmin);
        }

    }

    private void setupGeo() {

        this.setupGeoEC();
        this.setupGeoFTOF2();
        this.setupGeoFTOF1B();

    }

    private void setupGeoEC() {

        this.ECMatcher = new CrossMatcherToECClusters(this, this.idEC);
        this.ECMatcher.setupGeo();
        this.ECMatcher.setMinDistance(this.minDistance_ECAL);

    }

    private void setupGeoFTOF2() {

        this.FTOF2Matcher = new CrossMatcherToFTOFHits(this, 3);
        this.FTOF2Matcher.setupGeo();
        this.FTOF2Matcher.setMinDistanceCSI(this.minDistanceCSI_FTOF2);
        this.FTOF2Matcher.setMinDistanceY(this.minDistanceY_FTOF2);

    }

    private void setupGeoFTOF1B() {

        this.FTOF1BMatcher = new CrossMatcherToFTOFHits(this, 2);
        this.FTOF1BMatcher.setupGeo();
        this.FTOF1BMatcher.setMinDistanceCSI(this.minDistanceCSI_FTOF1B);
        this.FTOF1BMatcher.setMinDistanceY(this.minDistanceY_FTOF1B);

    }

    private void setupGeoFTOF1A() {

        this.FTOF1AMatcher = new CrossMatcherToFTOFHits(this, 1);
        this.FTOF1AMatcher.setupGeo();
        this.FTOF1AMatcher.setMinDistanceCSI(this.minDistanceCSI_FTOF1A);
        this.FTOF1AMatcher.setMinDistanceY(this.minDistanceY_FTOF1A);

    }

    private void setupHistograms() {

        this.allH1F = new ArrayList<H1F>();
        this.allH2F = new ArrayList<H2F>();
        h1_minCrossClusterDistance = new ArrayList<H1F>();

        for (int isector = 1; isector <= 6; isector++) {
            h1_minCrossClusterDistance.add(new H1F("h1_minCrossClusterDistance_" + isector, "minCrossClusterDistance_" + isector + ";distance - cm", 100, 0, 500));
            allH1F.add(h1_minCrossClusterDistance.get(isector - 1));
        }

        h1_allClustersE = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allClustersE.add(new H1F("h1_allClustersE_" + isector, "h1_allClustersE_" + isector + ";Energy (GeV)", 100, 0, .3));
            allH1F.add(h1_allClustersE.get(isector - 1));
        }

        h1_closerClustersE = new ArrayList<H1F>();
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_closerClustersE.add(new H1F("h1_closerClustersE_" + isector, "h1_closerClustersE_" + isector + ";Energy (GeV)", 100, 0, .3));
            allH1F.add(h1_closerClustersE.get(isector - 1));
        }

        h2_ThetaPhiAllGEN = new H2F("h2_ThetaPhiAllGEN", "h2_ThetaPhiAllGEN", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllGEN);
        h2_ThetaPhiAllREC = new H2F("h2_ThetaPhiAllREC", "h2_ThetaPhiAllREC", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllREC);
        h2_ThetaPhiAllTRIGGER = new H2F("h2_ThetaPhiAllTRIGGER", "h2_ThetaPhiAllTRIGGER", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllTRIGGER);
        h2_ThetaPhiAllTRIGGER2 = new H2F("h2_ThetaPhiAllTRIGGER2", "h2_ThetaPhiAllTRIGGER2", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllTRIGGER2);

        h2_ThetaPhiQPGEN = new H2F("h2_ThetaPhiQPGEN", "h2_ThetaPhiQPGEN", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPGEN);
        h2_ThetaPhiQPREC = new H2F("h2_ThetaPhiQPREC", "h2_ThetaPhiQPREC", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPGEN);
        h2_ThetaPhiQPTRIGGER = new H2F("h2_ThetaPhiQPTRIGGER", "h2_ThetaPhiQPTRIGGER", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPGEN);
        h2_ThetaPhiQPTRIGGER2 = new H2F("h2_ThetaPhiQPTRIGGER2", "h2_ThetaPhiQPTRIGGER2", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPGEN);

        h2_ThetaPhiQMGEN = new H2F("h2_ThetaPhiQMGEN", "h2_ThetaPhiQMGEN", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMGEN);
        h2_ThetaPhiQMREC = new H2F("h2_ThetaPhiQMREC", "h2_ThetaPhiQMREC", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMGEN);
        h2_ThetaPhiQMTRIGGER = new H2F("h2_ThetaPhiQMTRIGGER", "h2_ThetaPhiQMTRIGGER", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMGEN);
        h2_ThetaPhiQMTRIGGER2 = new H2F("h2_ThetaPhiQMTRIGGER2", "h2_ThetaPhiQMTRIGGER2", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMGEN);

        h1_vsDistanceAllREC = new H1F("h1_vsDistanceAllREC", "h1_vsDistanceAllREC", 400, 0., 400.);
        allH1F.add(h1_vsDistanceAllREC);
        h1_vsDistanceQPREC = new H1F("h1_vsDistanceQPREC", "h1_vsDistanceQPREC", 400, 0., 400.);
        allH1F.add(h1_vsDistanceQPREC);
        h1_vsDistanceQMREC = new H1F("h1_vsDistanceQMREC", "h1_vsDistanceQMREC", 400, 0., 400.);
        allH1F.add(h1_vsDistanceQMREC);

        h1_vsDistance1TRIGGER = new H1F("h1_vsDistance1TRIGGER", "h1_vsDistance3TRIGGER", 400, 0., 400.);
        allH1F.add(h1_vsDistance1TRIGGER);
        h1_vsDistance2TRIGGER = new H1F("h1_vsDistance2TRIGGER", "h1_vsDistance2TRIGGER", 400, 0., 400.);
        allH1F.add(h1_vsDistance2TRIGGER);
        h1_vsDistance3TRIGGER = new H1F("h1_vsDistance3TRIGGER", "h1_vsDistance3TRIGGER", 400, 0., 400.);
        allH1F.add(h1_vsDistance3TRIGGER);

        h1_vsDistance1TRIGGER2 = new H1F("h1_vsDistance1TRIGGER2", "h1_vsDistance3TRIGGER2", 400, 0., 400.);
        allH1F.add(h1_vsDistance1TRIGGER2);
        h1_vsDistance2TRIGGER2 = new H1F("h1_vsDistance2TRIGGER2", "h1_vsDistance2TRIGGER2", 400, 0., 400.);
        allH1F.add(h1_vsDistance2TRIGGER2);
        h1_vsDistance3TRIGGER2 = new H1F("h1_vsDistance3TRIGGER2", "h1_vsDistance3TRIGGER2", 400, 0., 400.);
        allH1F.add(h1_vsDistance3TRIGGER2);

        h1_vsDistanceAllREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllTRIGGER_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllTRIGGER2_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllEFF_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllEFF2_vsP = new ArrayList<H1F>();

        h1_vsDistanceQPREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPTRIGGER_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPTRIGGER2_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPEFF_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPEFF2_vsP = new ArrayList<H1F>();

        h1_vsDistanceQMREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMTRIGGER_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMTRIGGER2_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMEFF_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMEFF2_vsP = new ArrayList<H1F>();

        h2_ThetaPhiAllGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllTRIGGER_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllTRIGGER2_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllEFF_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllEFF2_vsP = new ArrayList<H2F>();

        h2_ThetaPhiQPGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPTRIGGER_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPTRIGGER2_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPEFF_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPEFF2_vsP = new ArrayList<H2F>();

        h2_ThetaPhiQMGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMTRIGGER_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMTRIGGER2_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMEFF_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMEFF2_vsP = new ArrayList<H2F>();

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {

            h1_vsDistanceAllREC_vsP.add(new H1F("h1_vsDistanceAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPREC_vsP.add(new H1F("h1_vsDistanceQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMREC_vsP.add(new H1F("h1_vsDistanceQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            h1_vsDistanceAllTRIGGER_vsP.add(new H1F("h1_vsDistanceAllTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceAllTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPTRIGGER_vsP.add(new H1F("h1_vsDistanceQPTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQPTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMTRIGGER_vsP.add(new H1F("h1_vsDistanceQMTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQMTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            h1_vsDistanceAllTRIGGER2_vsP.add(new H1F("h1_vsDistanceAllTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceAllTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPTRIGGER2_vsP.add(new H1F("h1_vsDistanceQPTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQPTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMTRIGGER2_vsP.add(new H1F("h1_vsDistanceQMTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQMTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            allH1F.add(h1_vsDistanceAllREC_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQPREC_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQMREC_vsP.get(ibin));
            allH1F.add(h1_vsDistanceAllTRIGGER_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQPTRIGGER_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQMTRIGGER_vsP.get(ibin));
            allH1F.add(h1_vsDistanceAllTRIGGER2_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQPTRIGGER2_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQMTRIGGER2_vsP.get(ibin));

            h2_ThetaPhiAllGEN_vsP.add(new H2F("h2_ThetaPhiAllGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiAllGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiAllREC_vsP.add(new H2F("h2_ThetaPhiAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiAllTRIGGER_vsP.add(new H2F("h2_ThetaPhiAllTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiAllTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiAllTRIGGER2_vsP.add(new H2F("h2_ThetaPhiAllTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiAllTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));

            allH2F.add(h2_ThetaPhiAllGEN_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiAllREC_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiAllTRIGGER_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiAllTRIGGER2_vsP.get(ibin));

            h2_ThetaPhiQPGEN_vsP.add(new H2F("h2_ThetaPhiQPGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQPGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQPREC_vsP.add(new H2F("h2_ThetaPhiQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQPTRIGGER_vsP.add(new H2F("h2_ThetaPhiQPTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQPTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQPTRIGGER2_vsP.add(new H2F("h2_ThetaPhiQPTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQPTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));

            allH2F.add(h2_ThetaPhiQPGEN_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQPREC_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQPTRIGGER_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQPTRIGGER2_vsP.get(ibin));

            h2_ThetaPhiQMGEN_vsP.add(new H2F("h2_ThetaPhiQMGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQMGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQMREC_vsP.add(new H2F("h2_ThetaPhiQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQMTRIGGER_vsP.add(new H2F("h2_ThetaPhiQMTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQMTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQMTRIGGER2_vsP.add(new H2F("h2_ThetaPhiQMTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQMTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));

            allH2F.add(h2_ThetaPhiQMGEN_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQMREC_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQMTRIGGER_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQMTRIGGER2_vsP.get(ibin));

        }

        /* MIP: ~ 2 MeV / cm, TOF2 thick=5 cm -> MPI=10 MeV */
        h2_FTOF2EnergyAll_LR = new H2F("h2_FTOF2EnergyAll_LR", "h2_FTOF2EnergyAll_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF2EnergyAll_LR);

        h2_FTOF2EnergyMatched_LR = new H2F("h2_FTOF2EnergyMatched_LR", "h2_FTOF2EnergyMatched_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF2EnergyMatched_LR);

        h2_FTOF1BEnergyAll_LR = new H2F("h2_FTOF1BEnergyAll_LR", "h2_FTOF1BEnergyAll_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1BEnergyAll_LR);

        h2_FTOF1BEnergyMatched_LR = new H2F("h2_FTOF1BEnergyMatched_LR", "h2_FTOF1BEnergyMatched_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1BEnergyMatched_LR);

        h2_FTOF1AEnergyAll_LR = new H2F("h2_FTOF1AEnergyAll_LR", "h2_FTOF1AEnergyAll_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1AEnergyAll_LR);

        h2_FTOF1AEnergyMatched_LR = new H2F("h2_FTOF1AEnergyMatched_LR", "h2_FTOF1AEnergyMatched_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1AEnergyMatched_LR);

    }

    private void processHistograms() {
        /* Last operations */
        h1_vsDistanceAllEFF = h1_vsDistanceAllREC.histClone("h1_vsDistanceAllEFF");
        h1_vsDistanceAllEFF.setTitle("h1_vsDistanceAllEFF");
        h1_vsDistanceAllEFF.divide(1. * nTracksWithR3Cross);

        h1_vsDistanceQPEFF = h1_vsDistanceQPREC.histClone("h1_vsDistanceQPEFF");
        h1_vsDistanceQPEFF.setTitle("h1_vsDistanceQPEFF");
        h1_vsDistanceQPEFF.divide(1. * nTracksQPWithR3Cross);

        h1_vsDistanceQMEFF = h1_vsDistanceQMREC.histClone("h1_vsDistanceQMEFF");
        h1_vsDistanceQMEFF.setTitle("h1_vsDistanceQMEFF");
        h1_vsDistanceQMEFF.divide(1. * nTracksQMWithR3Cross);

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {
            h1_vsDistanceAllEFF_vsP.add(h1_vsDistanceAllREC_vsP.get(ibin).histClone("h1_vsDistanceAllEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1)));
            h1_vsDistanceAllEFF_vsP.get(ibin).setTitle("h1_vsDistanceAllEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1));
            h1_vsDistanceAllEFF_vsP.get(ibin).divide(1. * nTracks_vsP.get(ibin));

            h1_vsDistanceQPEFF_vsP.add(h1_vsDistanceQPREC_vsP.get(ibin).histClone("h1_vsDistanceQPEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1)));
            h1_vsDistanceQPEFF_vsP.get(ibin).setTitle("h1_vsDistanceQPEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1));
            h1_vsDistanceQPEFF_vsP.get(ibin).divide(1. * nTracksQP_vsP.get(ibin));

            h1_vsDistanceQMEFF_vsP.add(h1_vsDistanceQMREC_vsP.get(ibin).histClone("h1_vsDistanceQMEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1)));
            h1_vsDistanceQMEFF_vsP.get(ibin).setTitle("h1_vsDistanceQMEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1));
            h1_vsDistanceQMEFF_vsP.get(ibin).divide(1. * nTracksQM_vsP.get(ibin));
        }

        h2_ThetaPhiAllEFF = H2F.divide(h2_ThetaPhiAllTRIGGER2, h2_ThetaPhiAllREC);
        h2_ThetaPhiQPEFF = H2F.divide(h2_ThetaPhiQPTRIGGER2, h2_ThetaPhiQPREC);
        h2_ThetaPhiQMEFF = H2F.divide(h2_ThetaPhiQMTRIGGER2, h2_ThetaPhiQMREC);

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {
            h2_ThetaPhiAllEFF_vsP.add(H2F.divide(h2_ThetaPhiAllTRIGGER2_vsP.get(ibin), h2_ThetaPhiAllREC_vsP.get(ibin)));
            h2_ThetaPhiQPEFF_vsP.add(H2F.divide(h2_ThetaPhiQPTRIGGER2_vsP.get(ibin), h2_ThetaPhiQPREC_vsP.get(ibin)));
            h2_ThetaPhiQMEFF_vsP.add(H2F.divide(h2_ThetaPhiQMTRIGGER2_vsP.get(ibin), h2_ThetaPhiQMREC_vsP.get(ibin)));

        }

        h1_vsDistance1EFF = h1_vsDistance1TRIGGER.histClone("h1_vsDistance1EFF");
        h1_vsDistance1EFF.setTitle("h1_vsDistance1EFF");
        h1_vsDistance1EFF.divide(1. * nEvents);

        h1_vsDistance2EFF = h1_vsDistance2TRIGGER.histClone("h1_vsDistance2EFF");
        h1_vsDistance2EFF.setTitle("h1_vsDistance2EFF");
        h1_vsDistance2EFF.divide(1. * nEvents);

        h1_vsDistance3EFF = h1_vsDistance3TRIGGER.histClone("h1_vsDistance3EFF");
        h1_vsDistance3EFF.setTitle("h1_vsDistance3EFF");
        h1_vsDistance3EFF.divide(1. * nEvents);

        h1_vsDistance1EFF2 = h1_vsDistance1TRIGGER2.histClone("h1_vsDistance1EFF2");
        h1_vsDistance1EFF2.setTitle("h1_vsDistance1EFF2");
        h1_vsDistance1EFF2.divide(1. * nEvents);

        h1_vsDistance2EFF2 = h1_vsDistance2TRIGGER2.histClone("h1_vsDistance2EFF2");
        h1_vsDistance2EFF2.setTitle("h1_vsDistance2EFF2");
        h1_vsDistance2EFF2.divide(1. * nEvents);

        h1_vsDistance3EFF2 = h1_vsDistance3TRIGGER2.histClone("h1_vsDistance3EFF2");
        h1_vsDistance3EFF2.setTitle("h1_vsDistance3EFF2");
        h1_vsDistance3EFF2.divide(1. * nEvents);

    }

    public H1F getHistogram1D(String name) {
        for (H1F h : this.allH1F) {
            if (name.equals(h.getName()) == true) {
                return h;
            }
        }
        System.out.println("histogram: " + name + " not found");
        return null;
    }

    public H2F getHistogram2D(String name) {
        for (H2F h : this.allH2F) {
            if (name.equals(h.getName()) == true) {
                return h;
            }
        }
        System.out.println("histogram: " + name + " not found");
        return null;
    }

    private void showHistograms() {

        TCanvas trk1 = new TCanvas("trk1", 1600, 1000);
        trk1.divide(2, 2);
        trk1.cd(1);
        for (int isector = 1; isector <= 6; isector++) {
            h1_minCrossClusterDistance.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_minCrossClusterDistance.get(isector - 1));
            } else {
                trk1.draw(h1_minCrossClusterDistance.get(isector - 1), "same");
            }
        }

        trk1.cd(2);
        for (int isector = 1; isector <= 6; isector++) {
            h1_allClustersE.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_allClustersE.get(isector - 1));
            } else {
                trk1.draw(h1_allClustersE.get(isector - 1), "same");
            }
        }

        trk1.cd(3);
        for (int isector = 1; isector <= 6; isector++) {
            h1_closerClustersE.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_closerClustersE.get(isector - 1));
            } else {
                trk1.draw(h1_closerClustersE.get(isector - 1), "same");
            }
        }

        TCanvas trkAllP = new TCanvas("allP", 1600, 1000);
        trkAllP.divide(5, 3);
        trkAllP.cd(0);
        trkAllP.draw(h1_vsDistanceAllEFF);
        h1_vsDistanceQPEFF.setLineColor(2);
        trkAllP.draw(h1_vsDistanceQPEFF, "same");
        h1_vsDistanceQMEFF.setLineColor(3);
        trkAllP.draw(h1_vsDistanceQMEFF, "same");

        trkAllP.cd(5);
        trkAllP.draw(h2_ThetaPhiQPGEN, "colz");

        trkAllP.cd(6);
        trkAllP.draw(h2_ThetaPhiQPREC, "colz");

        trkAllP.cd(7);
        trkAllP.draw(h2_ThetaPhiQPTRIGGER, "colz");

        trkAllP.cd(8);
        trkAllP.draw(h2_ThetaPhiQPTRIGGER2, "colz");

        trkAllP.cd(9);
        trkAllP.getCanvas().getPad(9).getAxisZ().setLog(true);
        trkAllP.draw(h2_ThetaPhiQPEFF, "colz");

        trkAllP.cd(10);
        trkAllP.draw(h2_ThetaPhiQMGEN, "colz");

        trkAllP.cd(11);
        trkAllP.draw(h2_ThetaPhiQMREC, "colz");

        trkAllP.cd(12);
        trkAllP.draw(h2_ThetaPhiQMTRIGGER, "colz");

        trkAllP.cd(13);
        trkAllP.draw(h2_ThetaPhiQMTRIGGER2, "colz");

        trkAllP.cd(14);
        trkAllP.getCanvas().getPad(14).getAxisZ().setLog(true);
        trkAllP.draw(h2_ThetaPhiQMEFF, "colz");

        ArrayList<TCanvas> TCanvasArray = new ArrayList<TCanvas>();
        for (int ii = 0; ii < Parray.size() - 1; ii++) {
            TCanvasArray.add(new TCanvas("P_" + Parray.get(ii) + "_" + Parray.get(ii + 1), 1600, 1000));
            TCanvasArray.get(ii).divide(5, 3);

            TCanvasArray.get(ii).cd(0);
            TCanvasArray.get(ii).draw(h1_vsDistanceAllEFF_vsP.get(ii));
            h1_vsDistanceQPEFF_vsP.get(ii).setLineColor(2);
            TCanvasArray.get(ii).draw(h1_vsDistanceQPEFF_vsP.get(ii), "same");
            h1_vsDistanceQMEFF_vsP.get(ii).setLineColor(3);
            TCanvasArray.get(ii).draw(h1_vsDistanceQMEFF_vsP.get(ii), "same");

            TCanvasArray.get(ii).cd(5);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPGEN_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(6);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPREC_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(7);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPTRIGGER_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(8);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPTRIGGER2_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(9);
            TCanvasArray.get(ii).getCanvas().getPad(9).getAxisZ().setLog(true);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPEFF_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(10);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMGEN_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(11);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMREC_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(12);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMTRIGGER_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(13);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMTRIGGER2_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(14);
            TCanvasArray.get(ii).getCanvas().getPad(14).getAxisZ().setLog(true);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMEFF_vsP.get(ii), "colz");

        }

        TCanvas ceff = new TCanvas("ceff", 1600, 1000);
        ceff.divide(2, 3);

        ceff.cd(0);
        ceff.draw(h1_vsDistance1EFF);
        h1_vsDistance2EFF.setLineColor(2);
        ceff.draw(h1_vsDistance2EFF, "same");
        h1_vsDistance3EFF.setLineColor(3);
        ceff.draw(h1_vsDistance3EFF, "same");

        ceff.cd(1);
        ceff.draw(h1_vsDistance1TRIGGER);
        h1_vsDistance2TRIGGER.setLineColor(2);
        ceff.draw(h1_vsDistance2TRIGGER, "same");
        h1_vsDistance3TRIGGER.setLineColor(3);
        ceff.draw(h1_vsDistance3TRIGGER, "same");

        ceff.cd(2);
        ceff.draw(h1_vsDistance1EFF2);
        h1_vsDistance2EFF2.setLineColor(2);
        ceff.draw(h1_vsDistance2EFF2, "same");
        h1_vsDistance3EFF2.setLineColor(3);
        ceff.draw(h1_vsDistance3EFF2, "same");

        ceff.cd(3);
        ceff.draw(h1_vsDistance1TRIGGER2);
        h1_vsDistance2TRIGGER2.setLineColor(2);
        ceff.draw(h1_vsDistance2TRIGGER2, "same");
        h1_vsDistance3TRIGGER2.setLineColor(3);
        ceff.draw(h1_vsDistance3TRIGGER2, "same");

        TCanvas cTOF2 = new TCanvas("cTOF2", 1600, 1600);
        cTOF2.divide(2, 3);
        cTOF2.cd(0);
        cTOF2.getCanvas().getPad(0).getAxisZ().setLog(true);
        cTOF2.draw(h2_FTOF2EnergyAll_LR, "colz");

        cTOF2.cd(1);
        cTOF2.getCanvas().getPad(1).getAxisZ().setLog(true);
        cTOF2.draw(h2_FTOF2EnergyMatched_LR, "colz");
        // cTOF2.draw(h2tmp,"colz");

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
            if (genParticlesDB == null) continue;
            nGenParticles = dataReader.makeGeneratedParticles(genParticlesDB, genParticles);

            /* Make tracks and crosses */
            DataBank tracksData = event.getBank("HitBasedTrkg::HBTracks");
            DataBank crossesData = event.getBank("HitBasedTrkg::HBCrosses");
            // DataBank tracksData = event.getBank("TimeBasedTrkg::TBTracks");
            // DataBank crossesData = event.getBank("TimeBasedTrkg::TBCrosses");
            if (tracksData == null) continue;
            if (crossesData == null) continue;

            tracks.clear();
            crosses.clear();

            dataReader.makeTriggerTracks(tracksData, crossesData, tracks, crosses);

            /*
             * Check how many MC particles have been reconstructed and match MC
             * - Track
             */
            int nReconstructedParticles = dataReader.matchReconstructedTracks(genParticles, tracks);

            if (tracksData.rows() >= 1) nEvents1Rec++;
            if (tracksData.rows() >= 2) nEvents2Rec++;
            if (tracksData.rows() >= 3) nEvents3Rec++;

            /* Get calorimeter clusters - later will search for PCAL */
            DataBank clustersEC = event.getBank("ECAL::clusters");
            if (clustersEC == null) continue;
            int nClustersEC = dataReader.makeCaloClusters(idEC, clustersEC, clusters);

            /*
             * Get FTOF hits - later will search for panel 2 (and maybe for
             * panel 1A)
             */
            DataBank hitsFTOF = event.getBank("FTOF::rawhits");
            if (hitsFTOF == null) continue;
            int nHitsFTOF2 = dataReader.readFTOFHits(hitsFTOF, hitsFTOF2, hitsFTOF1B, hitsFTOF1A);

            /* Perform matchings */
            /* Crosses */
            this.matchToClusters(crosses, clusters);
            this.matchToFTOFHits(3, crosses, hitsFTOF2);
            this.matchToFTOFHits(2, crosses, hitsFTOF1B);

            /* Tracks */
            this.matchToClusters(tracks, clusters);
            this.matchToFTOFHits(3, tracks, hitsFTOF2);
            this.matchToFTOFHits(2, tracks, hitsFTOF1B);
            
            
            /*Here do the "time-window analysis"
             * 
             */
            for (int ii=0;ii<AnalysisClass.nSectors_CLAS12;ii++){
                /*1A-1B coinc.*/
                this.matchFTOF1result[ii]=this.hitsMatcherFTOF1[ii].MatchHits(hitsFTOF1A[ii], hitsFTOF1B[ii]);             
            }
            /*
             * First, do a "direct" analysis - meaningfull if there's exactly
             * ONE generated particle
             */

            // if (nGenParticles == 1) {
            this.doDirectAnalysis();
            // }
            /*
             * Then, do an "inverse" analysis - basically what the trigger does:
             * take all R3 crosses - doesn't matter if matched to a track or
             * not, match them with ECAL.
             */

            this.doInverseCrossAnalysis();

        } /* End loop on events */

        this.processHistograms();
        if (this.doShowHistograms == true) {
            this.showHistograms();
        }

    }

}
