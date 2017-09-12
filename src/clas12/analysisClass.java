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

    /* Matching */
    private double minDistance = 50.;

    /* Momentum array */
    ArrayList<Double> Parray; // Index of lower-bound momenta (a part from
                              // latest, that is last bin
                              // higher-bound momentum)

    /* Variables */
    int nevent = -1;
    int nCrosses = 0;
    int nTracks = 0;
    int nTracksQP = 0;
    int nTracksQM = 0;

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
    H2F h2_ThetaPhiAllEFF;

    H2F h2_ThetaPhiQPGEN;
    H2F h2_ThetaPhiQPREC;
    H2F h2_ThetaPhiQPEFF;

    H2F h2_ThetaPhiQMGEN;
    H2F h2_ThetaPhiQMREC;
    H2F h2_ThetaPhiQMEFF;

    ArrayList<H2F> h2_ThetaPhiAllGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiAllREC_vsP;
    ArrayList<H2F> h2_ThetaPhiAllEFF_vsP;

    ArrayList<H2F> h2_ThetaPhiQPGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiQPREC_vsP;
    ArrayList<H2F> h2_ThetaPhiQPEFF_vsP;

    ArrayList<H2F> h2_ThetaPhiQMGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiQMREC_vsP;
    ArrayList<H2F> h2_ThetaPhiQMEFF_vsP;

    H1F h1_vsDistanceAllREC;
    H1F h1_vsDistanceAllEFF;

    H1F h1_vsDistanceQPREC;
    H1F h1_vsDistanceQPEFF;

    H1F h1_vsDistanceQMREC;
    H1F h1_vsDistanceQMEFF;

    ArrayList<H1F> h1_vsDistanceAllREC_vsP;
    ArrayList<H1F> h1_vsDistanceAllEFF_vsP;

    ArrayList<H1F> h1_vsDistanceQPREC_vsP;
    ArrayList<H1F> h1_vsDistanceQPEFF_vsP;

    ArrayList<H1F> h1_vsDistanceQMREC_vsP;
    ArrayList<H1F> h1_vsDistanceQMEFF_vsP;

    H1F h1_vsDistance1REC;
    H1F h1_vsDistance1EFF;
    H1F h1_vsDistance2REC;
    H1F h1_vsDistance2EFF;
    H1F h1_vsDistance3REC;
    H1F h1_vsDistance3EFF;

    H1F h1_vsDistance1MatchedCrossesREC;
    H1F h1_vsDistance1MatchedCrossesEFF;
    H1F h1_vsDistance2MatchedCrossesREC;
    H1F h1_vsDistance2MatchedCrossesEFF;
    H1F h1_vsDistance3MatchedCrossesREC;
    H1F h1_vsDistance3MatchedCrossesEFF;

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

    private int makeGeneratedParticles(DataBank genParticlesDB, List<Particle> genParticles) {
        int nGenParticles = 0;

        int nrows = genParticlesDB.rows();
        for (int loop = 0; loop < nrows; loop++) {
            Particle genParticle = new Particle(genParticlesDB.getInt("pid", loop), genParticlesDB.getFloat("px", loop), genParticlesDB.getFloat("py", loop),
                    genParticlesDB.getFloat("pz", loop), genParticlesDB.getFloat("vx", loop), genParticlesDB.getFloat("vy", loop), genParticlesDB.getFloat(
                            "vz", loop));
            genParticles.add(genParticle);
            nGenParticles++;
        }
        return nGenParticles;
    }

    private int makeCaloClusters(int detId, DataBank clustersDataBank, List<ECCluster>[] clusters) {
        int nClusters = clustersDataBank.rows();
        int nClustersDet = 0;

        for (int loop = 0; loop < nClusters; loop++) {
            int layer = clustersDataBank.getByte("layer", loop);
            if (layer / 3 != detId) continue; // 0: PCAL, 1:EC-IN, 2:EC-OUT

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

            if (energy < this.clusterEMin) continue; /*
                                                      * just consider clusters
                                                      * above my thr
                                                      */

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
     * (momentum of the track is in CLAS frame!)
     */
    private void makeTriggerTracks(DataBank tracksHBDataBank, DataBank crossesHBDataBank, List<TriggerTrack> tracks, List<Cross> crosses) {
        /* First, make all tracks */
        int nTracks = tracksHBDataBank.rows();
        int nCrosses = crossesHBDataBank.rows();

        for (int itrack = 0; itrack < nTracks; itrack++) {
            int sector = tracksHBDataBank.getByte("sector", itrack);
            int crossID3 = tracksHBDataBank.getShort("Cross3_ID", itrack);
            if (crossID3 == -1) continue;

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
            track.set_Id(analysisClass.crossIsAssociatedToTrack); /*
                                                                   * Use id to
                                                                   * save this
                                                                   */

            /* Add it to the list of crosses - TriggerTrack extends Cross */
            crosses.add(track);

        }
        /* Now make all the other R3 crosses */

        for (int loop = 0; loop < nCrosses; loop++) {
            int region = crossesHBDataBank.getByte("region", loop);
            if (region != 3) continue;
            short status = crossesHBDataBank.getShort("status", loop);


            /*
             * Already done before for crosses associated to tracks, need to
             * avoid double-counting
             */
            if (status == analysisClass.crossIsAssociatedToTrack) {
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
            cross.set_Id(0);

            /* Add elements to the lists */
            crosses.add(cross);

        }

    }

    public MatchingResult matchCrossToClusters(Cross cross, List<ECCluster> clusters) {

        MatchingResult result = new MatchingResult();
        result.setResult(0);
        result.setMinDistance(99999.);

        if (clusters.size() == 0) return result;

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
            System.out.println("The given cross does not intersect EC plane?" + pCross.toString() + " " + vCross.toString());
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
        if (dmin < this.minDistance) {
            this.h1_closerClustersE.get(sector - 1).fill(clusters.get(imin).energy);
        }
        h1_minCrossClusterDistance.get(sector - 1).fill(dmin);

        result.setResult(1); /* At least 1 cluster was found */
        result.setMinDistance(dmin);

        return result;
    }

    private void matchTracksToClusters(List<TriggerTrack> tracks, List<ECCluster> clusters[]) {

        ArrayList<MatchingResult> theMatchingResults = new ArrayList<MatchingResult>();

        for (TriggerTrack track : tracks) {

            double theta = Math.toDegrees(track.getMomentum().theta());
            double phi = Math.toDegrees(track.getMomentum().phi());
            double p = track.getMomentum().mag();
            int imom = this.getMomentumIndex(p);
            if (imom < 0) continue;
            byte q = track.getCharge();
            Integer value;

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

            nTracks++;

            value = nTracks_vsP.get(imom);
            value = value + 1; // increment value
            nTracks_vsP.set(imom, value); // replace value

            if (q > 0) {
                nTracksQP++;

                value = nTracksQP_vsP.get(imom);
                value = value + 1; // increment value
                nTracksQP_vsP.set(imom, value); // replace value

            }
            if (q < 0) {
                nTracksQM++;

                value = nTracksQM_vsP.get(imom);
                value = value + 1; // increment value
                nTracksQM_vsP.set(imom, value); // replace value

            }
            /*
             * Do the matching between this track and the clusters in this
             * sector
             */
            MatchingResult matching = matchCrossToClusters(track, clusters[track.get_Sector() - 1]);
            theMatchingResults.add(matching);

            // No ECclusters at all have been found in this
            // sector;
            if (matching.getResult() == 0)
                continue;
            else {
                double dMin = matching.getMinDistance();

                if (dMin < this.minDistance) {
                    h2_ThetaPhiAllREC.fill(phi, theta);
                    h2_ThetaPhiAllREC_vsP.get(imom).fill(phi, theta);

                    if (q > 0) {
                        h2_ThetaPhiQPREC.fill(phi, theta);
                        h2_ThetaPhiQPREC_vsP.get(imom).fill(phi, theta);
                    }
                    if (q < 0) {
                        h2_ThetaPhiQMREC.fill(phi, theta);
                        h2_ThetaPhiQMREC_vsP.get(imom).fill(phi, theta);
                    }

                }

                for (int ibin = 0; ibin < h1_vsDistanceAllREC.getxAxis().getNBins(); ibin++) {
                    double thisD = h1_vsDistanceAllREC.getxAxis().getBinCenter(ibin);
                    if (dMin < thisD) {
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
    }

    private void matchCrossesToClusters(List<Cross> crosses, List<ECCluster> clusters[]) {
        ArrayList<Double> matchingDistances = new ArrayList<Double>();
        ArrayList<Double> matchingDistancesMatchedCrosses = new ArrayList<Double>();
        double dMin, thisD;
        int nMatchings;
        for (Cross cross : crosses) {
            MatchingResult result = matchCrossToClusters(cross, clusters[cross.get_Sector() - 1]);
            nCrosses++;
            if (result.getResult() == 0)
                continue; // No ECclusters at all have been found in this
                          // sector;
            else {
                dMin = result.getMinDistance();
                matchingDistances.add(new Double(dMin));
                if (cross.get_Id() == analysisClass.crossIsAssociatedToTrack) {
                    matchingDistancesMatchedCrosses.add(new Double(dMin));
                }
            }
        }
        /*
         * Now, matchingDistances contains, for each cross in this event that
         * has been matched, the corresponding distance
         */
        int toPrint = 1;
        for (int ibin = 0; ibin < h1_vsDistance1REC.getxAxis().getNBins(); ibin++) {
            thisD = h1_vsDistance1REC.getxAxis().getBinCenter(ibin);
            nMatchings = 0;
            for (Double matching : matchingDistances) {
                if (matching < thisD) nMatchings++;
            }
            if (nMatchings >= 1) h1_vsDistance1REC.incrementBinContent(ibin);
            if (nMatchings >= 2) h1_vsDistance2REC.incrementBinContent(ibin);
            if (nMatchings >= 3) h1_vsDistance3REC.incrementBinContent(ibin);
        }

        /*
         * Above histograms were made considering ALL crosses. I repeat,
         * considering only crosses that are matched to a track
         */
        for (int ibin = 0; ibin < h1_vsDistance1REC.getxAxis().getNBins(); ibin++) {
            thisD = h1_vsDistance1REC.getxAxis().getBinCenter(ibin);
            nMatchings = 0;
            for (Double matching : matchingDistancesMatchedCrosses) {
                if (matching < thisD) nMatchings++;
            }
            if (nMatchings >= 1) h1_vsDistance1MatchedCrossesREC.incrementBinContent(ibin);
            if (nMatchings >= 2) h1_vsDistance2MatchedCrossesREC.incrementBinContent(ibin);
            if (nMatchings >= 3) h1_vsDistance3MatchedCrossesREC.incrementBinContent(ibin);
        }

    }

    /*
     * This method takes as input the MC bank - with the generated particles -
     * and the reconstructed tracks. For each generated particle, it tries to
     * figure out if that particle has been reconstructed properly. The method
     * returns the number of reconstructed tracks. To check if a particle has
     * been reconstructed:
     * 
     * 1) Loop over reconstructed particles in the same sector only --->NO. The
     * sector computed from gen. particle momentum may be different from the
     * actual particle sector due to solenoidal field!!! 2) Consider the same
     * charge 3) Check delta momentum 4) Check delta theta 5) Check delta phi
     * (this is a "sector-like" check, but the reconstructed track considers the
     * phi angle properly!)
     */
    private int checkReconstructedTracks(List<Particle> genParticles, List<TriggerTrack> recParticles) {
        int nReconstructed = 0;

        int sector, charge;
        double phi, P, theta;

        for (Particle genParticle : genParticles) {
            phi = genParticle.phi();
            P = genParticle.p();
            theta = genParticle.theta();
            charge = genParticle.charge();

            /*if (phi < 0) phi = phi + 2 * Math.PI;
            phi = phi + Math.toRadians(analysisClass.phiAngle / 2);
            if (phi > 2 * Math.PI) phi = phi - 2 * Math.PI;
            sector = (int) (Math.toDegrees(phi) / analysisClass.phiAngle);
            sector = sector + 1;*/

        

            for (TriggerTrack recParticle : recParticles) {

            
                /* Select same sector, same charge */
                // if (recParticle.get_Sector()!=sector) continue; //A.C.,
                // sector is computed from gen. momentum, and does not consider
                // sol. field
                
                if (recParticle.getCharge() != charge) continue;
                /* Very basic cuts over momentum and theta */
                if (Math.abs(recParticle.getMomentum().mag() - P) > 0.4) continue;
                if (Math.abs(recParticle.getMomentum().theta() - theta) > Math.toRadians(10.)) continue;
                if (Math.abs(recParticle.getMomentum().phi() - phi) > Math.toRadians(40.)) continue;
          
                
                
           
                recParticle.setGenParticle(genParticle);
                nReconstructed++;
                break; // if we arrive here - it means the matching was found.
                       // No reason to move forward in the loop

            }
        }

        return nReconstructed;
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

            if (nevent % 10000 == 0) System.out.println("Analyzed " + nevent + " events" + " " + nmultitrk);

            /* Get generated (MC) particles */
            DataBank genParticlesDB = event.getBank("MC::Particle");
            if (genParticlesDB == null) continue;
            List<Particle> genParticles = new ArrayList<Particle>();
            makeGeneratedParticles(genParticlesDB, genParticles);

            /* Get HB tracks */
            DataBank tracksHB = event.getBank("HitBasedTrkg::HBTracks");
            // DataBank tracksHB = event.getBank("TimeBasedTrkg::TBTracks");
            if (tracksHB == null) continue;

            /* Get HB crosses */
            DataBank crossesHB = event.getBank("HitBasedTrkg::HBCrosses");
            // DataBank crossesHB = event.getBank("TimeBasedTrkg::TBCrosses");
            if (crossesHB == null) continue;

            /* Make tracks and crosses */
            List<TriggerTrack> tracks = new ArrayList<TriggerTrack>();
            List<Cross> crosses = new ArrayList<Cross>();
            makeTriggerTracks(tracksHB, crossesHB, tracks, crosses);

            /* Check how many MC particles have been reconstructed */
            int nReconstructedParticles = checkReconstructedTracks(genParticles, tracks);
       

            if (nReconstructedParticles >= 1) nEvents1Rec++;
            if (nReconstructedParticles >= 2) nEvents2Rec++;
            if (nReconstructedParticles >= 3) nEvents3Rec++;

            /* Get calorimeter clusters - later will search for PCAL */
            DataBank clustersEC = event.getBank("ECAL::clusters");

            if (clustersEC == null) continue;

            List<ECCluster>[] clusters = new ArrayList[analysisClass.nSectors];
            for (int ii = 0; ii < analysisClass.nSectors; ii++) {
                clusters[ii] = new ArrayList<ECCluster>();
            }

            int nClustersEC = makeCaloClusters(idEC, clustersEC, clusters);

            /*
             * Now, I have the list of tracks and the list of clusters - sector
             * by sector. Need to match them
             */

            /*
             * First, do a "direct" analysis: select tracks, get their R3 cross,
             * and check if there's a matching with ECAL.This analysis is aimed
             * at maximizing the trigger efficiency for the signal
             */
            this.matchTracksToClusters(tracks, clusters);

            /*
             * Then, do an "inverse" analysis - basically what the trigger does:
             * take all R3 crosses - doesn't matter if matched to a track or
             * not, match them with ECAL.
             */
            this.matchCrossesToClusters(crosses, clusters);
        } /* End loop on events */

        this.processHistograms();
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
            } else if (splited[0].contains("minDistance")) {
                this.minDistance = Double.parseDouble(splited[1]);
            }
        }
        this.setupMomentum();
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
        Vector3D n = new Vector3D(Math.sin(Math.toRadians(analysisClass.thetaAngle)), 0., Math.cos(Math.toRadians(analysisClass.thetaAngle)));
        Point3D p = new Point3D(0, 0, this.l0EC / Math.cos(Math.toRadians(analysisClass.thetaAngle)));
        planeEC = new Plane3D(p, n);

    }

    private void setupHistograms() {
        h1_minCrossClusterDistance = new ArrayList<H1F>();

        for (int isector = 1; isector <= 6; isector++) {
            h1_minCrossClusterDistance.add(new H1F("h1_minCrossClusterDistance_" + isector, "minCrossClusterDistance_" + isector + ";distance - cm", 100, 0,
                    500));
        }

        h1_allClustersE = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allClustersE.add(new H1F("h1_allClustersE_" + isector, "h1_allClustersE_" + isector + ";Energy (GeV)", 100, 0, .3));
        }

        h1_closerClustersE = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_closerClustersE.add(new H1F("h1_cloaserClustersE_" + isector, "h1_closerClustersE_" + isector + ";Energy (GeV)", 100, 0, .3));
        }

        h2_ThetaPhiAllGEN = new H2F("h2_ThetaPhiAllGEN", "h2_ThetaPhiAllGEN", 100, -180., 180., 100, 0, 45.);
        h2_ThetaPhiAllREC = new H2F("h2_ThetaPhiAllREC", "h2_ThetaPhiAllREC", 100, -180., 180., 100, 0, 45.);

        h2_ThetaPhiQPGEN = new H2F("h2_ThetaPhiQPGEN", "h2_ThetaPhiQPGEN", 100, -180., 180., 100, 0, 45.);
        h2_ThetaPhiQPREC = new H2F("h2_ThetaPhiQPREC", "h2_ThetaPhiQPREC", 100, -180., 180., 100, 0, 45.);

        h2_ThetaPhiQMGEN = new H2F("h2_ThetaPhiQMGEN", "h2_ThetaPhiQMGEN", 100, -180., 180., 100, 0, 45.);
        h2_ThetaPhiQMREC = new H2F("h2_ThetaPhiQMREC", "h2_ThetaPhiQMREC", 100, -180., 180., 100, 0, 45.);
        // h2_ThetaPhiAllEFF = new H2F("h2_ThetaPhiAllEFF", "h2_ThetaPhiAllEFF",
        // 100,
        // -180., 180., 100, 0, 45.);

        h1_vsDistanceAllREC = new H1F("h1_vsDistanceAllREC", "h1_vsDistanceAllREC", 400, 0., 400.);
        h1_vsDistanceQPREC = new H1F("h1_vsDistanceQPREC", "h1_vsDistanceQPREC", 400, 0., 400.);
        h1_vsDistanceQMREC = new H1F("h1_vsDistanceQMREC", "h1_vsDistanceQMREC", 400, 0., 400.);

        h1_vsDistance1REC = new H1F("h1_vsDistance1REC", "h1_vsDistance3REC", 400, 0., 400.);
        h1_vsDistance2REC = new H1F("h1_vsDistance2REC", "h1_vsDistance2REC", 400, 0., 400.);
        h1_vsDistance3REC = new H1F("h1_vsDistance3REC", "h1_vsDistance3REC", 400, 0., 400.);

        h1_vsDistance1MatchedCrossesREC = new H1F("h1_vsDistance1MatchedCrossesREC", "h1_vsDistance1MatchedCrossesREC", 400, 0., 400.);
        h1_vsDistance2MatchedCrossesREC = new H1F("h1_vsDistance2MatchedCrossesREC", "h1_vsDistance2MatchedCrossesREC", 400, 0., 400.);
        h1_vsDistance3MatchedCrossesREC = new H1F("h1_vsDistance3MatchedCrossesREC", "h1_vsDistance3MatchedCrossesREC", 400, 0., 400.);

        h1_vsDistanceAllREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllEFF_vsP = new ArrayList<H1F>();

        h1_vsDistanceQPREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPEFF_vsP = new ArrayList<H1F>();

        h1_vsDistanceQMREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMEFF_vsP = new ArrayList<H1F>();

        h2_ThetaPhiAllGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllEFF_vsP = new ArrayList<H2F>();

        h2_ThetaPhiQPGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPEFF_vsP = new ArrayList<H2F>();

        h2_ThetaPhiQMGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMEFF_vsP = new ArrayList<H2F>();

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {
            h1_vsDistanceAllREC_vsP.add(new H1F("h1_vsDistanceAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceAllEFF_"
                    + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPREC_vsP.add(new H1F("h1_vsDistanceQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQPEFF_" + Parray.get(ibin)
                    + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMREC_vsP.add(new H1F("h1_vsDistanceQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h1_vsDistanceQMEFF_" + Parray.get(ibin)
                    + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            h2_ThetaPhiAllGEN_vsP.add(new H2F("h2_ThetaPhiAllGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiAllGEN_" + Parray.get(ibin)
                    + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 45.));
            h2_ThetaPhiAllREC_vsP.add(new H2F("h2_ThetaPhiAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiAllREC_" + Parray.get(ibin)
                    + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 45.));

            h2_ThetaPhiQPGEN_vsP.add(new H2F("h2_ThetaPhiQPGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQPGEN_" + Parray.get(ibin) + "_"
                    + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 45.));
            h2_ThetaPhiQPREC_vsP.add(new H2F("h2_ThetaPhiQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQPREC_" + Parray.get(ibin) + "_"
                    + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 45.));

            h2_ThetaPhiQMGEN_vsP.add(new H2F("h2_ThetaPhiQMGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQMGEN_" + Parray.get(ibin) + "_"
                    + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 45.));
            h2_ThetaPhiQMREC_vsP.add(new H2F("h2_ThetaPhiQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), "h2_ThetaPhiQMREC_" + Parray.get(ibin) + "_"
                    + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 45.));

        }

    }

    private void processHistograms() {
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

        h2_ThetaPhiAllEFF = H2F.divide(h2_ThetaPhiAllREC, h2_ThetaPhiAllGEN);
        h2_ThetaPhiQPEFF = H2F.divide(h2_ThetaPhiQPREC, h2_ThetaPhiQPGEN);
        h2_ThetaPhiQMEFF = H2F.divide(h2_ThetaPhiQMREC, h2_ThetaPhiQMGEN);

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {
            h2_ThetaPhiAllEFF_vsP.add(H2F.divide(h2_ThetaPhiAllREC_vsP.get(ibin), h2_ThetaPhiAllGEN_vsP.get(ibin)));
            h2_ThetaPhiQPEFF_vsP.add(H2F.divide(h2_ThetaPhiQPREC_vsP.get(ibin), h2_ThetaPhiQPGEN_vsP.get(ibin)));
            h2_ThetaPhiQMEFF_vsP.add(H2F.divide(h2_ThetaPhiQMREC_vsP.get(ibin), h2_ThetaPhiQMGEN_vsP.get(ibin)));
        }

        h1_vsDistance1EFF = h1_vsDistance1REC.histClone("h1_vsDistance1EFF");
        h1_vsDistance1EFF.setTitle("h1_vsDistance1EFF");
        h1_vsDistance1EFF.divide(1. * nEvents);

        h1_vsDistance2EFF = h1_vsDistance2REC.histClone("h1_vsDistance2EFF");
        h1_vsDistance2EFF.setTitle("h1_vsDistance2EFF");
        h1_vsDistance2EFF.divide(1. * nEvents);

        h1_vsDistance3EFF = h1_vsDistance3REC.histClone("h1_vsDistance3EFF");
        h1_vsDistance3EFF.setTitle("h1_vsDistance3EFF");
        h1_vsDistance3EFF.divide(1. * nEvents);

        h1_vsDistance1MatchedCrossesEFF = h1_vsDistance1MatchedCrossesREC.histClone("h1_vsDistance1MatchedCrossesEFF");
        h1_vsDistance1MatchedCrossesEFF.setTitle("h1_vsDistance1MatchedCrossesEFF");
        h1_vsDistance1MatchedCrossesEFF.divide(1. * nEvents1Rec);

        h1_vsDistance2MatchedCrossesEFF = h1_vsDistance2MatchedCrossesREC.histClone("h1_vsDistance2MatchedCrossesEFF");
        h1_vsDistance2MatchedCrossesEFF.setTitle("h1_vsDistance2MatchedCrossesEFF");
        h1_vsDistance2MatchedCrossesEFF.divide(1. * nEvents2Rec);

        h1_vsDistance3MatchedCrossesEFF = h1_vsDistance3MatchedCrossesREC.histClone("h1_vsDistance3MatchedCrossesEFF");
        h1_vsDistance3MatchedCrossesEFF.setTitle("h1_vsDistance3MatchedCrossesEFF");
        h1_vsDistance3MatchedCrossesEFF.divide(1. * nEvents3Rec);

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
        trkAllP.divide(3, 3);
        trkAllP.cd(0);
        trkAllP.draw(h1_vsDistanceAllEFF);
        h1_vsDistanceQPEFF.setLineColor(2);
        trkAllP.draw(h1_vsDistanceQPEFF, "same");
        h1_vsDistanceQMEFF.setLineColor(3);
        trkAllP.draw(h1_vsDistanceQMEFF, "same");

        trkAllP.cd(3);
        trkAllP.draw(h2_ThetaPhiQPGEN, "colz");

        trkAllP.cd(4);
        trkAllP.draw(h2_ThetaPhiQPREC, "colz");

        trkAllP.cd(5);
        trkAllP.draw(h2_ThetaPhiQPEFF, "colz");

        trkAllP.cd(6);
        trkAllP.draw(h2_ThetaPhiQMGEN, "colz");

        trkAllP.cd(7);
        trkAllP.draw(h2_ThetaPhiQMREC, "colz");

        trkAllP.cd(8);
        trkAllP.draw(h2_ThetaPhiQMEFF, "colz");

        ArrayList<TCanvas> TCanvasArray = new ArrayList<TCanvas>();
        for (int ii = 0; ii < Parray.size() - 1; ii++) {
            TCanvasArray.add(new TCanvas("P_" + Parray.get(ii) + "_" + Parray.get(ii + 1), 1600, 1000));
            TCanvasArray.get(ii).divide(3, 3);

            TCanvasArray.get(ii).cd(0);
            TCanvasArray.get(ii).draw(h1_vsDistanceAllEFF_vsP.get(ii));
            h1_vsDistanceQPEFF_vsP.get(ii).setLineColor(2);
            TCanvasArray.get(ii).draw(h1_vsDistanceQPEFF_vsP.get(ii), "same");
            h1_vsDistanceQMEFF_vsP.get(ii).setLineColor(3);
            TCanvasArray.get(ii).draw(h1_vsDistanceQMEFF_vsP.get(ii), "same");

            TCanvasArray.get(ii).cd(3);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPGEN_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(4);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPREC_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(5);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPEFF_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(6);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMGEN_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(7);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMREC_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(8);
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
        ceff.draw(h1_vsDistance1REC);
        h1_vsDistance2REC.setLineColor(2);
        ceff.draw(h1_vsDistance2REC, "same");
        h1_vsDistance3REC.setLineColor(3);
        ceff.draw(h1_vsDistance3REC, "same");

        ceff.cd(2);
        ceff.draw(h1_vsDistance1MatchedCrossesEFF);
        h1_vsDistance2MatchedCrossesEFF.setLineColor(2);
        ceff.draw(h1_vsDistance2MatchedCrossesEFF, "same");
        h1_vsDistance3MatchedCrossesEFF.setLineColor(3);
        ceff.draw(h1_vsDistance3MatchedCrossesEFF, "same");

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
            nTracks_vsP.add(new Integer(0));
            nTracksQP_vsP.add(new Integer(0));
            nTracksQM_vsP.add(new Integer(0));
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

}
