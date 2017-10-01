package clas12;

import java.util.List;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;

public class DataReader {

	private static final int ID_FTOF = 17;
	private static final int ID_CAL = 16;

	private AnalysisClass analysisClass;

	private double minClusterE_ECAL = 0.01; // GeV
	private double minE_FTOF2 = 0.2; // MeV
	private double minE_FTOF1B = 0.2; // MeV
	private double minE_FTOF1A = 0.2; // MeV

	private double minE_FTOF; // used only for rate histogramming

	boolean hasDCsegmentsL1[];
	boolean hasDCsegmentsL2[];

	public DataReader(AnalysisClass ana) {
		analysisClass = ana;
		hasDCsegmentsL1 = new boolean[AnalysisClass.nSectors_CLAS12];
		hasDCsegmentsL2 = new boolean[AnalysisClass.nSectors_CLAS12];
	}

	public void setMinE_FTOF(double minE_FTOF) {
		this.minE_FTOF = minE_FTOF;
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

	public int makeGeneratedParticles(DataBank genParticlesDB, List<MatchedParticle> genParticles) {
		int nGenParticles = 0;
		int nrows = genParticlesDB.rows();

		genParticles.clear();
		for (int loop = 0; loop < nrows; loop++) {
			MatchedParticle genParticle = new MatchedParticle(genParticlesDB.getInt("pid", loop), genParticlesDB.getFloat("px", loop),
					genParticlesDB.getFloat("py", loop), genParticlesDB.getFloat("pz", loop), genParticlesDB.getFloat("vx", loop),
					genParticlesDB.getFloat("vy", loop), genParticlesDB.getFloat("vz", loop));
			genParticles.add(genParticle);
			nGenParticles++;
		}
		return nGenParticles;
	}

	public int makeReconstructedParticles(DataBank recParticlesDB, DataBank recParticlesTOFDB, DataBank recParticlesCALDB, List<MatchedParticle> recParticles) {
		int nRecParticles = 0;
		int nrows = recParticlesDB.rows();
		int nrowsTOF = recParticlesTOFDB.rows();
		int nrowsCAL = recParticlesCALDB.rows();

		int detector, particle;
		boolean hasTOF, hasCAL;

		recParticles.clear();
		for (int loop = 0; loop < nrows; loop++) { /*
													 * Check if there's a corresponding entry in recParticlesTOFDB or in
													 * recParticlesCALDB
													 */

			hasTOF = false;
			hasCAL = false;

			for (int loopTOF = 0; loopTOF < nrowsTOF; loopTOF++) { /* Start with FTOF */
				detector = recParticlesTOFDB.getByte("detector", loopTOF);
				if (detector != DataReader.ID_FTOF) continue;
				particle = recParticlesTOFDB.getShort("pindex", loopTOF);
				if (particle == loop) {
					hasTOF = true;
					break;
				}
			}
			if (!hasTOF) {
				for (int loopCAL = 0; loopCAL < nrowsCAL; loopCAL++) {
					detector = recParticlesCALDB.getByte("detector", loopCAL);
					if (detector != DataReader.ID_CAL) continue;
					particle = recParticlesCALDB.getShort("pindex", loopCAL);
					if (particle == loop) {
						hasCAL = true;
						break;
					}
				}
			}

			if (hasTOF || hasCAL) {
				/*
				 * so far, PID is not good in recon. However, I will work with: -> Momentum
				 * (from tracking) -> Charge (from tracking) Hence, if q>0, I set pip (211), if
				 * q<0, I set pim (-211), if q=0, I set gamma (22)
				 */
				byte charge = recParticlesDB.getByte("charge", loop);
				int pid;
				if (charge > 0)
					pid = 211;
				else if (charge < 0)
					pid = -211;
				else
					pid = 22;
				MatchedParticle recParticle = new MatchedParticle(pid, 0.1, charge, recParticlesDB.getFloat("px", loop), recParticlesDB.getFloat("py", loop),
						recParticlesDB.getFloat("pz", loop), recParticlesDB.getFloat("vx", loop), recParticlesDB.getFloat("vy", loop),
						recParticlesDB.getFloat("vz", loop));
				recParticles.add(recParticle);

				nRecParticles++;
			}
		}
		return nRecParticles;
	}

	public int makeCaloClusters(int detId, DataBank clustersDataBank, List<ECCluster>[] clusters, int DC, boolean[] hasDCsegments) {
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

			analysisClass.getHistogram1D("h1_allPCAL_E_" + sector).fill(energy);

			if (energy < this.minClusterE_ECAL) continue; /*
															 * just consider clusters above my thr
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

			/*
			 * Add it to the list of clusters - if there's any DC requirement, then don't
			 * put these in the output array!
			 */
			if ((DC == 0) || (hasDCsegments[sector - 1])) {
				clusters[sector - 1].add(cluster);
			}
		}

		return nClustersDet;
	}

	public int makeFTOFHits(DataBank hitsRawFTOFDataBank, DataBank hitsReconFTOFDataBank, List<ReconTOFHit>[] hitsFTOF2, List<ReconTOFHit>[] hitsFTOF1B,
			List<ReconTOFHit>[] hitsFTOF1A, int DC, boolean hasDCsegments[]) {

		int nRawFTOFHits = hitsRawFTOFDataBank.rows();
		int nReconFTOFHits = hitsReconFTOFDataBank.rows();

		if (nRawFTOFHits != nReconFTOFHits) {
			System.out.println(
					"Warning! This event: " + analysisClass.nevent + "has different number of TOF hits Raw/Recon: " + nRawFTOFHits + " " + nReconFTOFHits);
		}

		for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
			hitsFTOF2[ii].clear();
			hitsFTOF1B[ii].clear();
			hitsFTOF1A[ii].clear();
		}

		for (int loop = 0; loop < nRawFTOFHits; loop++) {
			int sector = hitsReconFTOFDataBank.getByte("sector", loop);
			int layer = hitsReconFTOFDataBank.getByte("layer", loop);
			int component = hitsReconFTOFDataBank.getShort("component", loop);
			short id = hitsReconFTOFDataBank.getShort("id", loop);

			float energy = hitsReconFTOFDataBank.getFloat("energy", loop);
			float time = hitsReconFTOFDataBank.getFloat("time", loop);

			float x = hitsReconFTOFDataBank.getFloat("x", loop);
			float y = hitsReconFTOFDataBank.getFloat("y", loop);
			float z = hitsReconFTOFDataBank.getFloat("z", loop);

			/* Following lines are used to rotate the point to the sector system */
			Point3D p0 = new Point3D(x, y, z);
			p0.rotateZ(-(Math.toRadians((sector - 1) * AnalysisClass.phiAngle_CLAS12)));

			float energyL = hitsRawFTOFDataBank.getFloat("energy_left", loop);
			float energyR = hitsRawFTOFDataBank.getFloat("energy_right", loop);
			float timeL = hitsRawFTOFDataBank.getFloat("time_left", loop);
			float timeR = hitsRawFTOFDataBank.getFloat("time_right", loop);

			/* Should I put here also the L/R time coincidence? */
			switch (layer) {
			case 1:
				if ((energyL > this.minE_FTOF1A) && (energyR > this.minE_FTOF1A)) {

					ReconTOFHit hit = new ReconTOFHit(sector, layer, component, id, energyL, energyR, timeL, timeR, energy, time);
					hit.setP0(p0);

					if ((DC == 0) || (hasDCsegments[sector - 1])) {
						hitsFTOF1A[sector - 1].add(hit);
					}
					if (energy > this.minE_FTOF) {
						analysisClass.getHistogram1D("h1_rateFTOF1A_" + sector).fill(component);
					}
					analysisClass.getHistogram1D("h1_allTOF1A_E_" + sector).fill(energy);
					analysisClass.getHistogram2D("h2_FTOF1AEnergyAll_LR").fill(energyL, energyR);
				}
				break;
			case 2:
				if ((energyL > this.minE_FTOF1B) && (energyR > this.minE_FTOF1B)) {
					ReconTOFHit hit = new ReconTOFHit(sector, layer, component, id, energyL, energyR, timeL, timeR, energy, time);
					hit.setP0(p0);
					if ((DC == 0) || (hasDCsegments[sector - 1])) {
						hitsFTOF1B[sector - 1].add(hit);
					}
					if (energy > this.minE_FTOF) {
						analysisClass.getHistogram1D("h1_rateFTOF1B_" + sector).fill(component);
					}
					analysisClass.getHistogram1D("h1_allTOF1B_E_" + sector).fill(energy);
					analysisClass.getHistogram2D("h2_FTOF1BEnergyAll_LR").fill(energyL, energyR);
				}
				break;
			case 3: // panel2
				if ((energyL > this.minE_FTOF2) && (energyR > this.minE_FTOF2)) {
					ReconTOFHit hit = new ReconTOFHit(sector, layer, component, id, energyL, energyR, timeL, timeR, energy, time);
					hit.setP0(p0);
					if ((DC == 0) || (hasDCsegments[sector - 1])) {
						hitsFTOF2[sector - 1].add(hit);
					}
					if (energy > this.minE_FTOF) {
						analysisClass.getHistogram1D("h1_rateFTOF2_" + sector).fill(component);
					}
					analysisClass.getHistogram1D("h1_allTOF2_E_" + sector).fill(energy);
					analysisClass.getHistogram2D("h2_FTOF2EnergyAll_LR").fill(energyL, energyR);
				}
				break;
			}
		}
		return nReconFTOFHits;
	}

	/*
	 * Make trigger tracks and all-R3 tracks, in SECTOR reference frame (momentum of
	 * the track is in CLAS frame!)
	 */
	public void makeDCCrossesAndTracks(DataBank tracksHBDataBank, DataBank crossesHBDataBank, List<TrackMatchedToGen> tracks, List<MatchedCross> crosses) {
		/* First, make all tracks */
		int nTracks = tracksHBDataBank.rows();
		int nCrosses = crossesHBDataBank.rows();

		for (int itrack = 0; itrack < nTracks; itrack++) {
			int sector = tracksHBDataBank.getByte("sector", itrack);
			int crossID3 = tracksHBDataBank.getShort("Cross3_ID", itrack);
			if (crossID3 == -1) continue; /*
											 * Case when this track has no R3-cross associated with: I simply ignore this
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
			Vector3D u = new Vector3D(ux, uy, uz);
			Point3D p0 = new Point3D(x0, y0, z0);
			Vector3D eu = new Vector3D(eux, euy, euz);
			Point3D ep0 = new Point3D(ex0, ey0, ez0);

			p0.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
			u.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

			ep0.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));
			eu.rotateY(Math.toRadians(AnalysisClass.thetaAngleDC_CLAS12));

			/*
			 * Uncomment following lines to go to CLAS12 system
			 * p0.rotateZ(Math.toRadians((sector - 1) * phiAngle));
			 * p1.rotateZ(Math.toRadians((sector - 1) * phiAngle));
			 */
			track.set_Dir(u.toPoint3D());
			track.set_DirErr(eu.toPoint3D());
			track.set_Point(p0);
			track.set_PointErr(ep0);

			tracks.add(track);

			/*
			 * Mark the corresponding cross as associated to a track using the status
			 */
			crossesHBDataBank.setShort("status", crossID, AnalysisClass.crossIsAssociatedToTrack);
			track.set_Id(AnalysisClass.crossIsAssociatedToTrack); /*
																	 * Use id to save this
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
			 * Already done before for crosses associated to tracks, need to avoid
			 * double-counting
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

	public void makeDCSegments(DataBank segmentsDCDataBank, List<SimpleDCSegment> segments, boolean[] sectorHasR3Segments) {
		int nSegments = segmentsDCDataBank.rows();

		for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
			hasDCsegmentsL1[ii] = false;
			hasDCsegmentsL2[ii] = false;
			sectorHasR3Segments[ii] = false;
		}

		for (int isegment = 0; isegment < nSegments; isegment++) {
			int sector = segmentsDCDataBank.getByte("sector", isegment);
			int superlayer = segmentsDCDataBank.getByte("superlayer", isegment);

			segments.add(new SimpleDCSegment(sector, superlayer));

			if (superlayer == 5) hasDCsegmentsL1[sector - 1] = true;
			if (superlayer == 6) hasDCsegmentsL2[sector - 1] = true;

		}

		for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
			if (hasDCsegmentsL1[ii] && hasDCsegmentsL2[ii]) sectorHasR3Segments[ii] = true;
		}
	}

	/*
	 * This method takes as input the MC bank - with the generated particles - and
	 * the reconstructed tracks. For each generated particle, it tries to figure out
	 * if that particle has been reconstructed properly. The method returns the
	 * number of reconstructed tracks. To check if a particle has been
	 * reconstructed:
	 * 
	 * 1) Loop over reconstructed particles in the same sector only --->NO. The
	 * sector computed from gen. particle momentum may be different from the actual
	 * particle sector due to solenoidal field!!! 2) Consider the same charge 3)
	 * Check delta momentum 4) Check delta theta 5) Check delta phi (this is a
	 * "sector-like" check, but the reconstructed track considers the phi angle
	 * properly!)
	 */
	public int matchReconstructedTracks(List<MatchedParticle> genParticles, List<TrackMatchedToGen> recParticles) {
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
			 * Math.toRadians(analysisClass.phiAngle / 2); if (phi > 2 * Math.PI) phi = phi
			 * - 2 * Math.PI; sector = (int) (Math.toDegrees(phi) / analysisClass.phiAngle);
			 * sector = sector + 1;
			 */

			for (TrackMatchedToGen recParticle : recParticles) {

				/* Select same sector, same charge */
				// if (recParticle.get_Sector()!=sector) continue; //A.C.,
				// sector is computed from gen. momentum, and does not consider
				// sol. field

				/*
				 * CLAS12 DC specifics are deltaP/P < 1% deltaTheta < 1 mrad -> 0.057 deg
				 * deltaPhi < 1 mrad / sinTheta
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

	/* Check if generated particles have been reconstructed, if so match them */
	public void matchReconstructedParticles(List<MatchedParticle> recParticles, List<MatchedParticle> genParticles) {

		for (int iGen = 0; iGen < genParticles.size(); iGen++) {
			MatchedParticle genParticle = genParticles.get(iGen);

			for (int iRec = 0; iRec < recParticles.size(); iRec++) {
				MatchedParticle recParticle = recParticles.get(iRec);

				if (analysisClass.nevent == 17969) {
					System.out.println("GEN: " + genParticle.vector().vect().toString() + " REC: " + recParticle.vector().vect().toString());
				}

				if (recParticle.charge() != genParticle.charge()) continue;
				if (Math.abs(recParticle.p() - genParticle.p()) / genParticle.p() > 0.2) continue;
				if (Math.abs(recParticle.theta() - genParticle.theta()) > Math.toRadians(10.)) continue;
				if (Math.abs(recParticle.phi() - genParticle.phi()) > Math.toRadians(20.)) continue;

				genParticle.setMatched(true);
				genParticle.setIdMatch(iRec);

				recParticle.setMatched(true);
				recParticle.setIdMatch(iGen);
				break; // if we arrive here - it means the matching was found.
			}

		}

	}
}
