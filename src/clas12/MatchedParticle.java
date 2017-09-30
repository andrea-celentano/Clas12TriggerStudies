package clas12;

import org.jlab.clas.physics.Particle;

public class MatchedParticle extends Particle {

    private boolean isMatched;
    private int idMatch;

    public MatchedParticle() {
        super();
        isMatched = false;
        idMatch = -1;

    }

    public MatchedParticle(MatchedParticle p) {
        super(p);
        isMatched = false;
        idMatch = -1;
    }

    public MatchedParticle(int pid, double px, double py, double pz, double vx, double vy, double vz) {
        super(pid, px, py, pz, vx, vy, vz);
        isMatched = false;
        idMatch = -1;
    }

    public MatchedParticle(int pid, double px, double py, double pz) {
        super(pid, px, py, pz);
        isMatched = false;
        idMatch = -1;
    }

    public MatchedParticle(int pid, double mass, byte charge, double px, double py, double pz, double vx, double vy, double vz) {
        super(pid, mass, charge, px, py, pz, vx, vy, vz);
        isMatched = false;
        idMatch = -1;

    }

    public boolean isMatched() {
        return isMatched;
    }

    public void setMatched(boolean isMatched) {
        this.isMatched = isMatched;
    }

    public int getIdMatch() {
        return idMatch;
    }

    public void setIdMatch(int idMatch) {
        this.idMatch = idMatch;
    }
}
