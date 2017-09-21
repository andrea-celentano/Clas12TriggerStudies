package clas12;

import java.util.List;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;

public class DataReaderAndMatcher {

    private AnalysisClass analysisClass;

    private double minClusterE_ECAL = 0.01; // GeV
    private double minE_FTOF2 = 0.02; // MeV
    private double minE_FTOF1B = 0.02; // MeV
    private double minE_FTOF1A = 0.02; // MeV

    public DataReaderAndMatcher(AnalysisClass ana) {
        analysisClass = ana;
    }

    public double getMinClusterE_ECAL() {
        return minClusterE_ECAL;
    }

    public void setMinClusterE_ECAL(double minClusterE_ECAL) {
        this.minClusterE_ECAL = minClusterE_ECAL;
    }

    public double getMinE_FTOF1A() {
        return minE_FTOF1B;
    }

    public void setMinE_FTOF1A(double minE_FTOF1A) {
        this.minE_FTOF1A = minE_FTOF1A;
    }

    public double getMinE_FTOF1B() {
        return minE_FTOF1B;
    }

    public void setMinE_FTOF1B(double minE_FTOF1B) {
        this.minE_FTOF1B = minE_FTOF1B;
    }

    public double getMinE_FTOF2() {
        return minE_FTOF2;
    }

    public void setMinE_FTOF2(double minE_FTOF2) {
        this.minE_FTOF2 = minE_FTOF2;
    }

    public int makeGeneratedParticles(DataBank genParticlesDB, List<Particle> genParticles) {
        int nGenParticles = 0;
        int nrows = genParticlesDB.rows();

        genParticles.clear();
        for (int loop = 0; loop < nrows; loop++) {
            Particle genParticle = new Particle(genParticlesDB.getInt("pid", loop), genParticlesDB.getFloat("px", loop), genParticlesDB.getFloat("py", loop), genParticlesDB.getFloat("pz", loop), genParticlesDB.getFloat("vx", loop), genParticlesDB.getFloat("vy", loop), genParticlesDB.getFloat("vz", loop));
            genParticles.add(genParticle);
            nGenParticles++;
        }
        return nGenParticles;
    }

    public int makeCaloClusters(int detId, DataBank clustersDataBank, List<ECCluster>[] clusters) {
        int nClusters = clustersDataBank.rows();
        int nClustersDet = 0;

        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            clusters[ii].clear();
        }

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

            analysisClass.getHistogram1D("h1_allClustersE_" + sector).fill(energy);

            if (energy < this.minClusterE_ECAL) continue; /*
                                                           * just consider
                                                           * clusters above my
                                                           * thr
                                                           */

            /* Following lines are used to rotate the point to the sector system */
            Point3D p0 = new Point3D(x, y, z);
            p0.rotateZ(-(Math.toRadians((sector - 1) * AnalysisClass.phiAngle_CLAS12)));

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

    public int readFTOFHits(DataBank hitsFTOFDataBank, List<SimpleTOFHit>[] hitsFTOF2, List<SimpleTOFHit>[] hitsFTOF1B, List<SimpleTOFHit>[] hitsFTOF1A) {

        int nFTOFHits = hitsFTOFDataBank.rows();

        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            hitsFTOF2[ii].clear();
            hitsFTOF1B[ii].clear();
            hitsFTOF1A[ii].clear();
        }

        for (int loop = 0; loop < nFTOFHits; loop++) {
            int sector = hitsFTOFDataBank.getByte("sector", loop);
            int layer = hitsFTOFDataBank.getByte("layer", loop);
            int component = hitsFTOFDataBank.getShort("component", loop);
            short id = hitsFTOFDataBank.getShort("id", loop);
            float energyL = hitsFTOFDataBank.getFloat("energy_left", loop);
            float energyR = hitsFTOFDataBank.getFloat("energy_right", loop);
            float timeL = hitsFTOFDataBank.getFloat("time_left", loop);
            float timeR = hitsFTOFDataBank.getFloat("time_right", loop);

            SimpleTOFHit hit = new SimpleTOFHit(sector, layer, component, id, energyL, energyR, timeL, timeR);

            // System.out.println(sector+" "+layer+" "+component+" "+energyL+" "+energyR+" "+timeL+" "+timeR);

            switch (layer) {
            case 1:
                if ((energyL > this.minE_FTOF1A) && (energyR > this.minE_FTOF1A)) {
                    hitsFTOF1A[sector - 1].add(hit);
                    analysisClass.getHistogram2D("h2_FTOF1AEnergyAll_LR").fill(energyL, energyR);
                }
                break;
            case 2:
                if ((energyL > this.minE_FTOF1B) && (energyR > this.minE_FTOF1B)) {
                    hitsFTOF1B[sector - 1].add(hit);
                    analysisClass.getHistogram2D("h2_FTOF1BEnergyAll_LR").fill(energyL, energyR);
                }
                break;
            case 3: // panel2
                if ((energyL > this.minE_FTOF2) && (energyR > this.minE_FTOF2)) {
                    hitsFTOF2[sector - 1].add(hit);
                    analysisClass.getHistogram2D("h2_FTOF2EnergyAll_LR").fill(energyL, energyR);
                }
                break;
            }
        }
        return nFTOFHits;
    }

    /*
     * Make trigger tracks and all-R3 tracks, in SECTOR reference frame
     * (momentum of the track is in CLAS frame!)
     */
    public void makeTriggerTracks(DataBank tracksHBDataBank, DataBank crossesHBDataBank, List<TrackMatchedToGen> tracks, List<MatchedCross> crosses) {
        /* First, make all tracks */
        int nTracks = tracksHBDataBank.rows();
        int nCrosses = crossesHBDataBank.rows();

        for (int itrack = 0; itrack < nTracks; itrack++) {
            int sector = tracksHBDataBank.getByte("sector", itrack);
            int crossID3 = tracksHBDataBank.getShort("Cross3_ID", itrack);
            if (crossID3 == -1) continue; /*
                                           * Case when this track has no
                                           * R3-cross associated with: I simply
                                           * ignore this
                                           */

            byte charge = tracksHBDataBank.getByte("q", itrack);
            float px = tracksHBDataBank.getFloat("p0_x", itrack);
            float py = tracksHBDataBank.getFloat("p0_y", itrack);
            float pz = tracksHBDataBank.getFloat("p0_z", itrack);
            Vector3 p = new Vector3(px, py, pz);

            /* Make trigger track */
            TrackMatchedToGen track = new TrackMatchedToGen(sector);
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
                continue; // ignore this cross
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

            p0.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
            u.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

            ep0.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
            eu.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

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
            crossesHBDataBank.setShort("status", crossID, AnalysisClass.crossIsAssociatedToTrack);
            track.set_Id(AnalysisClass.crossIsAssociatedToTrack); /*
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
            if (status == AnalysisClass.crossIsAssociatedToTrack) {
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

            p0.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
            u.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

            ep0.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
            eu.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

            /*
             * These lines are commented because I work in the SECTOR system
             * p0.rotateZ(Math.toRadians((sector - 1) * phiAngle));
             * p1.rotateZ(Math.toRadians((sector - 1) * phiAngle));
             */

            /* Create the cross */
            MatchedCross cross = new MatchedCross(sector);
            cross.set_Dir(u.toPoint3D());
            cross.set_DirErr(eu.toPoint3D());
            cross.set_Point(p0);
            cross.set_PointErr(ep0);
            cross.set_Id(0);

            /* Add elements to the lists */
            crosses.add(cross);
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
    public int matchReconstructedTracks(List<Particle> genParticles, List<TrackMatchedToGen> recParticles) {
        int nMatched = 0;
        int iGen = 0;
        int charge;
        double phi, P, theta;

        for (Particle genParticle : genParticles) {
            phi = genParticle.phi();
            P = genParticle.p();
            theta = genParticle.theta();
            charge = genParticle.charge();

            /*
             * if (phi < 0) phi = phi + 2 * Math.PI; phi = phi +
             * Math.toRadians(analysisClass.phiAngle / 2); if (phi > 2 *
             * Math.PI) phi = phi - 2 * Math.PI; sector = (int)
             * (Math.toDegrees(phi) / analysisClass.phiAngle); sector = sector +
             * 1;
             */

            for (TrackMatchedToGen recParticle : recParticles) {

                /* Select same sector, same charge */
                // if (recParticle.get_Sector()!=sector) continue; //A.C.,
                // sector is computed from gen. momentum, and does not consider
                // sol. field

                /*
                 * CLAS12 DC specifics are deltaP/P < 1% deltaTheta < 1 mrad ->
                 * 0.057 deg deltaPhi < 1 mrad / sinTheta
                 * 
                 * -> I take much bigger values
                 */

                if (recParticle.getCharge() != charge) continue;
                /* Very basic cuts over momentum and theta */
                if (Math.abs(recParticle.getMomentum().mag() - P) / P > 0.3) continue;
                if (Math.abs(recParticle.getMomentum().theta() - theta) > Math.toRadians(10.)) continue;
                if (Math.abs(recParticle.getMomentum().phi() - phi) > Math.toRadians(40.)) continue;

                recParticle.setGenParticle(iGen, genParticle);
                nMatched++;
                break; // if we arrive here - it means the matching was found.
                       // No reason to move forward in the loop of reconstructed
            }

            iGen++;
        }

        return nMatched;
    }
}
