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
//import org.jlab.service.ec.ECCluster;

import java.util.List;
import java.util.ArrayList;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

public class analysisClass {

    /* Here goes common definitions */
    private static final Vector3D xAxis = new Vector3D(1., 0., 0.);
    private static final Vector3D yAxis = new Vector3D(0., 1., 0.);
    private static final Vector3D zAxis = new Vector3D(0., 0., 1.);

    private static double phiAngle = 60.;
    private static double thetaAngle = 25.;
    private static int nSectors = 6;

    /* EC */
    private static double l0PCAL = 697.7;
    private static double hPCAL = 385.1;
    private static double bPCAL = 394.;
    private static double dPCAL = 290.8;

    private static double l0ECIN;
    private static double l0ECOUT;

    private double l0EC; /*
                          * The distance between the target and the EC plane
                          * along the 25 deg line
                          */
    private double bEC, hEC; /* The EC base and height */
    private double dEC; /*
                         * The distance between the triangle vertex and the
                         * 25-deg line intersection point
                         */
    private Plane3D planeEC;

    /* Crosses */
    private static short crossIsAssociatedToTrack = 255;

    private String inputFileName = "";
    private int nEvents = -1;

    /* Clusters */
    private int idEC = 0; // 0: PCAL, 1: EC-IN, 2: EC-OUT
    private double clusterEMin = 0.;

    /* Here goes histograms */
    ArrayList<H1F> h1_minCrossClusterDistance;
    ArrayList<H1F> h1_allClustersE;
    ArrayList<H1F> h1_closerClustersE;

    H2F h2_ThetaPhiAllGEN;
    H2F h2_ThetaPhiAllREC;
    H2F h2_ThetaPhiAllEFF;

    H2F h2_ThetaPhiQPGEN;
    H2F h2_ThetaPhiQPREC;
    H2F h2_ThetaPhiQPEFF;

    H2F h2_ThetaPhiQMGEN;
    H2F h2_ThetaPhiQMREF;
    H2F h2_ThetaPhiQMEFF;

    H1F h1_vsDistanceAllREC;
    H1F h1_vsDistanceAllEFF;

    H1F h1_vsDistanceQPREC;
    H1F h1_vsDistanceQPEFF;

    H1F h1_vsDistanceQMREC;
    H1F h1_vsDistanceQMEFF;

    H1F h1_crossVsDistanceAllREC;
    H1F h1_crossVsDistanceAllEFF;

    /* my cluster helper class */

    class ECCluster {

        public int sector, layer;
        public double energy, time;
        public Point3D p0;

        public ECCluster() {
            p0 = new Point3D();
        }

        public double getEnergy() {
            return energy;
        }

    }

    public static void main(String[] args) throws IOException {
        // TODO Auto-generated method stub

        analysisClass myClass = new analysisClass();

        myClass.setup(args);

        myClass.print();
        myClass.run();

    }

    private int makeCaloClusters(int detId, DataBank clustersDataBank, List<ECCluster>[] clusters) {
        int nClusters = clustersDataBank.rows();
        int nClustersDet = 0;

        for (int loop = 0; loop < nClusters; loop++) {
            int layer = clustersDataBank.getByte("layer", loop);
            if (layer / 3 != detId)
                continue; // 0: PCAL, 1:EC-IN, 2:EC-OUT

            int sector = clustersDataBank.getByte("sector", loop);
            if ((sector < 1) || (sector > 6)) {
                System.out.println("Bad sector" + sector);
                continue;
            }
            double x = clustersDataBank.getFloat("x", loop);
            double y = clustersDataBank.getFloat("y", loop);
            double z = clustersDataBank.getFloat("z", loop);
            double energy = clustersDataBank.getFloat("energy", loop);
            double time = clustersDataBank.getFloat("time", loop);

            this.h1_allClustersE.get(sector - 1).fill(energy);

            if (energy < this.clusterEMin)
                continue;

            /* Following lines are used to rotate the point to the sector system */
            Point3D p0 = new Point3D(x, y, z);
            p0.rotateZ(-(Math.toRadians((sector - 1) * phiAngle)));

            /* Create the cluster */
            ECCluster cluster = new ECCluster();
            cluster.sector = sector;
            cluster.layer = layer;
            cluster.time = time;
            cluster.energy = energy;
            cluster.p0 = p0;

            /* Add it to the list of clusters */
            clusters[sector - 1].add(cluster);
        }

        return nClustersDet;
    }

    /*
     * Make trigger tracks and all-R3 tracks, in SECTOR reference frame
     * (momentum of the track is in CLAS frame!
     */
    private void makeTriggerTracks(DataBank tracksHBDataBank, DataBank crossesHBDataBank, List<TriggerTrack> tracks,
            List<Cross> crosses) {
        /* First, make all tracks */
        int nTracks = tracksHBDataBank.rows();
        int nCrosses = crossesHBDataBank.rows();

        for (int itrack = 0; itrack < nTracks; itrack++) {
            int sector = tracksHBDataBank.getByte("sector", itrack);
            int crossID3 = tracksHBDataBank.getShort("Cross3_ID", itrack);
            if (crossID3 == -1)
                continue;

            byte charge = tracksHBDataBank.getByte("q", itrack);
            float px = tracksHBDataBank.getFloat("p0_x", itrack);
            float py = tracksHBDataBank.getFloat("p0_y", itrack);
            float pz = tracksHBDataBank.getFloat("p0_z", itrack);
            Vector3 p = new Vector3(px, py, pz);

            /* Make trigger track */
            TriggerTrack track = new TriggerTrack(sector);
            track.setMomentum(p);
            track.setCharge(charge);

            /* Set cross variables */
            int crossID = -1;
            for (int icross = 0; icross < nCrosses; icross++) {
                int id = crossesHBDataBank.getShort("id", icross);
                if (id == crossID3) {
                    crossID = icross;
                    break;
                }
            }
            if (crossID == -1) {
                System.out.println("Error with cross indexing: " + crossID3 + " " + crossID);
            }

            double x0 = crossesHBDataBank.getFloat("x", crossID);
            double y0 = crossesHBDataBank.getFloat("y", crossID);
            double z0 = crossesHBDataBank.getFloat("z", crossID);

            double ux = crossesHBDataBank.getFloat("ux", crossID);
            double uy = crossesHBDataBank.getFloat("uy", crossID);
            double uz = crossesHBDataBank.getFloat("uz", crossID);

            double ex0 = crossesHBDataBank.getFloat("err_x", crossID);
            double ey0 = crossesHBDataBank.getFloat("err_y", crossID);
            double ez0 = crossesHBDataBank.getFloat("err_z", crossID);

            double eux = crossesHBDataBank.getFloat("err_ux", crossID);
            double euy = crossesHBDataBank.getFloat("err_uy", crossID);
            double euz = crossesHBDataBank.getFloat("err_uz", crossID);

            /* Translate point from tilted to sector coordinates */
            /* For direction vector, simply imagine it defines a second point. */
            Vector3D u = new Vector3D(ux, uy, uz);
            Point3D p0 = new Point3D(x0, y0, z0);
            Vector3D eu = new Vector3D(eux, euy, euz);
            Point3D ep0 = new Point3D(ex0, ey0, ez0);

            p0.rotateY(Math.toRadians(thetaAngle));
            u.rotateY(Math.toRadians(thetaAngle));

            ep0.rotateY(Math.toRadians(thetaAngle));
            eu.rotateY(Math.toRadians(thetaAngle));

            /*
             * These lines are commented because I work in the SECTOR system
             * p0.rotateZ(Math.toRadians((sector - 1) * phiAngle));
             * p1.rotateZ(Math.toRadians((sector - 1) * phiAngle));
             */
            track.set_Dir(u.toPoint3D());
            track.set_DirErr(eu.toPoint3D());
            track.set_Point(p0);
            track.set_PointErr(ep0);

            tracks.add(track);
            /*
             * Mark the corresponding cross as associated to a track using the
             * status
             */
            crossesHBDataBank.setShort("status", crossID, analysisClass.crossIsAssociatedToTrack);

            /* Add it to the list of crosses - TriggerTrack extends Cross */
            crosses.add(track);

        }
        /* Now make all the other R3 crosses */

        for (int loop = 0; loop < nCrosses; loop++) {
            int region = crossesHBDataBank.getByte("region", loop);
            if (region != 3)
                continue;
            short status = crossesHBDataBank.getShort("status", loop);
            if (status == analysisClass.crossIsAssociatedToTrack) { /*Already done before*/
                continue;
            }
            int sector = crossesHBDataBank.getByte("sector", loop);

            int id = crossesHBDataBank.getShort("id", loop);
            double x0 = crossesHBDataBank.getFloat("x", loop);
            double y0 = crossesHBDataBank.getFloat("y", loop);
            double z0 = crossesHBDataBank.getFloat("z", loop);

            double ux = crossesHBDataBank.getFloat("ux", loop);
            double uy = crossesHBDataBank.getFloat("uy", loop);
            double uz = crossesHBDataBank.getFloat("uz", loop);

            double ex0 = crossesHBDataBank.getFloat("err_x", loop);
            double ey0 = crossesHBDataBank.getFloat("err_y", loop);
            double ez0 = crossesHBDataBank.getFloat("err_z", loop);

            double eux = crossesHBDataBank.getFloat("err_ux", loop);
            double euy = crossesHBDataBank.getFloat("err_uy", loop);
            double euz = crossesHBDataBank.getFloat("err_uz", loop);
            /* Translate point from tilted to sector coordinates */
            /* For direction vector, simply imagine it defines a second point. */
            Vector3D u = new Vector3D(ux, uy, uz);
            Point3D p0 = new Point3D(x0, y0, z0);
            Vector3D eu = new Vector3D(eux, euy, euz);
            Point3D ep0 = new Point3D(ex0, ey0, ez0);

            p0.rotateY(Math.toRadians(thetaAngle));
            u.rotateY(Math.toRadians(thetaAngle));

            ep0.rotateY(Math.toRadians(thetaAngle));
            eu.rotateY(Math.toRadians(thetaAngle));

            /*
             * These lines are commented because I work in the SECTOR system
             * p0.rotateZ(Math.toRadians((sector - 1) * phiAngle));
             * p1.rotateZ(Math.toRadians((sector - 1) * phiAngle));
             */

            /* Create the cross */
            Cross cross = new Cross(sector, region, id);
            cross.set_Dir(u.toPoint3D());
            cross.set_DirErr(eu.toPoint3D());
            cross.set_Point(p0);
            cross.set_PointErr(ep0);

            /* Add elements to the lists */
            crosses.add(cross);

        }

    }

    public TriggerResult matchCrossToClusters(Cross cross, List<ECCluster> clusters) {

        
        
        
        TriggerResult result = new TriggerResult();
        result.setResult(0);
        result.setMinDistance(99999.);
        
        if (clusters.size() == 0)
            return result;

        
        /* Do the matching between this cross and the ECCLusters in this sector */
        int sector = cross.get_Sector();
  
        /*
         * Create a line passing by the cross point and with direction equal to
         * the cross direction
         */
        Point3D pCross = cross.get_Point(); /* passage point */
        Vector3D vCross = cross.get_Dir().toVector3D();/* Direction */
        Line3D rayCross = new Line3D(pCross, vCross);

        /* Determine the intersection between this line and the EC-plane */
        Point3D intersectCross = new Point3D();
        if ((planeEC.intersectionRay(rayCross, intersectCross) != 1)) {
            System.out.println("The given cross does not intersect EC plane?" + pCross.toString() + " "
                    + vCross.toString());
            return result;
        }

        /*
         * System.out.println("This cross is: "+cross.get_Point().toString
         * ()+" "+cross.get_Dir().toString());
         */

        /*
         * now, find the clusters in this sector, and compute the distance
         */
        ArrayList<Double> distances = new ArrayList<Double>();
        for (ECCluster cluster : clusters) {
            Point3D pCluster = cluster.p0;
            double d = pCluster.distance(intersectCross);
            distances.add(new Double(d));
        }

        /* Find the minimum distance */
        double dmin = 9999.;
        int imin = 0;
        for (Double d : distances) {
            if (d < dmin) {
                dmin = d;
                imin = distances.indexOf(d);
            }
        }
        if (dmin<50.){
            this.h1_closerClustersE.get(sector - 1).fill(clusters.get(imin).energy);
        }
        h1_minCrossClusterDistance.get(sector - 1).fill(dmin);
        
        result.setResult(1); /* At least 1 cluster was found */
        result.setMinDistance(dmin);

        return result;
    }

    public void run() {

        HipoDataSource reader = new HipoDataSource();
        reader.open(this.inputFileName);
        System.out.println("There are: " + reader.getSize() + " events in the file");
        H1F hi = new H1F("hi_ftof_energy_PA", "hi_ftof_energy_PA", 100, 0.0, 1.0);

        int nevent = -1;
        int nmultitrk = 0;

        if (nEvents == -1) {
            nEvents = reader.getSize();
            System.out.println("Reading ALL ");
        } else {
            System.out.println("Reading: " + nEvents);
        }

        int nCrosses = 0;
        int nTracks = 0;
        int nTracksQP = 0;
        int nTracksQM = 0;

        while (reader.hasEvent() == true && nevent < nEvents) {
            DataEvent event = reader.getNextEvent();
            nevent++;
            if (nevent % 10000 == 0)
                System.out.println("Analyzed " + nevent + " events" + " " + nmultitrk);
            // event.show();

            /* Get HB tracks */
            DataBank tracksHB = event.getBank("HitBasedTrkg::HBTracks");
            if (tracksHB == null)
                continue;

            /* Get HB crosses */
            DataBank crossesHB = event.getBank("HitBasedTrkg::HBCrosses");
            if (crossesHB == null)
                continue;

            /* Make tracks and crosses */
            List<TriggerTrack> tracks = new ArrayList<TriggerTrack>();
            List<Cross> crosses = new ArrayList<Cross>();
            makeTriggerTracks(tracksHB, crossesHB, tracks, crosses);

            /* Search for R3 crosses, and fill the corresponding vectors */
            /* This method returns coordinates in the SECTOR system */

            /* Get calorimeter clusters - later will search for PCAL */
            DataBank clustersEC = event.getBank("ECAL::clusters");

            if (clustersEC == null)
                continue;

            List<ECCluster>[] clusters = new ArrayList[analysisClass.nSectors];
            for (int ii = 0; ii < analysisClass.nSectors; ii++) {
                clusters[ii] = new ArrayList<ECCluster>();
            }

            int nClustersEC = makeCaloClusters(idEC, clustersEC, clusters);

            /*
             * Now, I have the list of tracks and the list of clusters - sector
             * by sector. Need to match them
             */
            for (TriggerTrack track : tracks) {

                double theta = Math.toDegrees(track.getMomentum().theta());
                double phi = Math.toDegrees(track.getMomentum().phi());
                double p = track.getMomentum().mag();
                byte q = track.getCharge();

                TriggerResult trigger = matchCrossToClusters(track, clusters[track.get_Sector() - 1]);

                /* Fill the "generated histograms" */
                h2_ThetaPhiAllGEN.fill(phi, theta);
                nTracks++;
                if (q > 0)
                    nTracksQP++;
                if (q < 0)
                    nTracksQM++;
                if (trigger.getResult() == 0)
                    continue; // No ECclusters at all have been found in this
                              // sector;
                else {

                    h2_ThetaPhiAllREC.fill(phi, theta);

                    double dMin = trigger.getMinDistance();
                    for (int ibin = 0; ibin < h1_vsDistanceAllREC.getxAxis().getNBins(); ibin++) {
                        double thisD = h1_vsDistanceAllREC.getxAxis().getBinCenter(ibin);
                        if (dMin < thisD) {
                            if (q > 0) {
                                h1_vsDistanceQPREC.incrementBinContent(ibin);
                            }
                            if (q < 0) {
                                h1_vsDistanceQMREC.incrementBinContent(ibin);
                            }
                            h1_vsDistanceAllREC.incrementBinContent(ibin);
                        }
                    }
                }
            }

            /*
             * Need to repeat the SAME operation using ALL crosses - that are
             * the objects actually seen by the trigger!
             */
            for (Cross cross : crosses) {
                TriggerResult trigger = matchCrossToClusters(cross, clusters[cross.get_Sector() - 1]);
                nCrosses++;
                if (trigger.getResult() == 0)
                    continue; // No ECclusters at all have been found in this
                              // sector;
                else {
                    double dMin = trigger.getMinDistance();
                    for (int ibin = 0; ibin < h1_crossVsDistanceAllREC.getxAxis().getNBins(); ibin++) {
                        double thisD = h1_crossVsDistanceAllREC.getxAxis().getBinCenter(ibin);
                        if (dMin < thisD) {
                            h1_crossVsDistanceAllREC.incrementBinContent(ibin);
                        }
                    }
                }
            }
        }/* End loop on events */

        /* Last operations */
        h1_vsDistanceAllEFF = h1_vsDistanceAllREC.histClone("h1_vsDistanceAllEFF");
        h1_vsDistanceAllEFF.setTitle("h1_vsDistanceAllEFF");
        h1_vsDistanceAllEFF.divide(1. * nTracks);

        h1_vsDistanceQPEFF = h1_vsDistanceQPREC.histClone("h1_vsDistanceQPEFF");
        h1_vsDistanceQPEFF.setTitle("h1_vsDistanceQPEFF");
        h1_vsDistanceQPEFF.divide(1. * nTracksQP);

        h1_vsDistanceQMEFF = h1_vsDistanceQMREC.histClone("h1_vsDistanceQMEFF");
        h1_vsDistanceQMEFF.setTitle("h1_vsDistanceQMEFF");
        h1_vsDistanceQMEFF.divide(1. * nTracksQM);

        h2_ThetaPhiAllEFF = H2F.divide(h2_ThetaPhiAllREC, h2_ThetaPhiAllGEN);

        
        h1_crossVsDistanceAllEFF = h1_crossVsDistanceAllREC.histClone("h1_crossVsDistanceAllEFF");
        h1_crossVsDistanceAllEFF.setTitle("h1_crossVsDistanceAllEFF");
        h1_crossVsDistanceAllEFF.divide(1. * nCrosses);

        
        
        
        this.showHistograms();

    }

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
            } else if (splited[0].contains("idEC")) {
                this.idEC = Integer.parseInt(splited[1]);
            } else if (splited[0].contains("clusterE")) {
                this.clusterEMin = Double.parseDouble(splited[1]);
            }
        }

        this.setupHistograms();
        this.setupGeo();
    }

    private void setupGeo() {
        switch (this.idEC) {
        case 0:
            this.l0EC = analysisClass.l0PCAL;
            this.bEC = analysisClass.bPCAL;
            this.hEC = analysisClass.hPCAL;
            this.dEC = analysisClass.dPCAL;
            break;
        case 1:
            this.l0EC = analysisClass.l0ECIN;
            break;
        case 2:
            this.l0EC = analysisClass.l0ECOUT;
            break;
        default:
            this.l0EC = analysisClass.l0PCAL;
            break;
        }

        /*
         * Define the EC-plane, in the sector system (i.e. the plane of sector
         * n.1)
         */
        Vector3D n = new Vector3D(Math.sin(Math.toRadians(analysisClass.thetaAngle)), 0., Math.cos(Math
                .toRadians(analysisClass.thetaAngle)));
        Point3D p = new Point3D(0, 0, this.l0EC / Math.cos(Math.toRadians(analysisClass.thetaAngle)));
        planeEC = new Plane3D(p, n);

    }

    private void setupHistograms() {
        h1_minCrossClusterDistance = new ArrayList<H1F>();

        for (int isector = 1; isector <= 6; isector++) {
            h1_minCrossClusterDistance.add(new H1F("h1_minCrossClusterDistance_" + isector, "minCrossClusterDistance_"
                    + isector + ";distance - cm", 100, 0, 500));
        }

        h1_allClustersE = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allClustersE.add(new H1F("h1_allClustersE_" + isector, "h1_allClustersE_" + isector + ";Energy (GeV)",
                    100, 0, .3));
        }

        h1_closerClustersE = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_closerClustersE.add(new H1F("h1_cloaserClustersE_" + isector, "h1_closerClustersE_" + isector
                    + ";Energy (GeV)", 100, 0, .3));
        }

        h2_ThetaPhiAllGEN = new H2F("h2_ThetaPhiAllGEN", "h2_ThetaPhiAllGEN", 100, -180., 180., 100, 0, 45.);
        h2_ThetaPhiAllREC = new H2F("h2_ThetaPhiAllREC", "h2_ThetaPhiAllREC", 100, -180., 180., 100, 0, 45.);
        h2_ThetaPhiAllEFF = new H2F("h2_ThetaPhiAllEFF", "h2_ThetaPhiAllEFF", 100, -180., 180., 100, 0, 45.);

        h1_vsDistanceAllREC = new H1F("h1_vsDistanceAllREC", "h1_vsDistanceAllEFF", 400, 0., 400.);
        h1_vsDistanceAllEFF = new H1F("h1_vsDistanceAllEFF", "h1_vsDistanceAllEFF", 400, 0., 400.);

        h1_vsDistanceQPREC = new H1F("h1_vsDistanceQPREC", "h1_vsDistanceQPEFF", 400, 0., 400.);
        h1_vsDistanceQPEFF = new H1F("h1_vsDistanceQPEFF", "h1_vsDistanceQPEFF", 400, 0., 400.);

        h1_vsDistanceQMREC = new H1F("h1_vsDistanceQMREC", "h1_vsDistanceQMEFF", 400, 0., 400.);
        h1_vsDistanceQMEFF = new H1F("h1_vsDistanceQMEFF", "h1_vsDistanceQMEFF", 400, 0., 400.);

        h1_crossVsDistanceAllREC = new H1F("h1_crossVsDistanceAllREC", "h1_crossVsDistanceAllEFF", 400, 0., 400.);
        h1_crossVsDistanceAllEFF = new H1F("h1_crossVsDistanceAllEFF", "h1_crossVsDistanceAllEFF", 400, 0., 400.);

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

        TCanvas trk2 = new TCanvas("trk2", 1600, 1000);
        trk2.divide(3, 2);
        trk2.cd(0);
        trk2.draw(h1_vsDistanceAllEFF);
        h1_vsDistanceQPEFF.setLineColor(2);
        trk2.draw(h1_vsDistanceQPEFF, "same");
        h1_vsDistanceQMEFF.setLineColor(3);
        trk2.draw(h1_vsDistanceQMEFF, "same");

        trk2.cd(3);
        trk2.draw(h2_ThetaPhiAllGEN, "colz");

        trk2.cd(4);
        trk2.draw(h2_ThetaPhiAllREC, "colz");

        trk2.cd(5);
        trk2.draw(h2_ThetaPhiAllEFF, "colz");

        TCanvas trk3 = new TCanvas("trk3", 1600, 1000);
        trk3.divide(2, 2);
        trk3.cd(0);
        trk3.draw(h1_crossVsDistanceAllEFF);
        
    }

}
