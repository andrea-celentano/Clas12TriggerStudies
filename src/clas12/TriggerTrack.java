package clas12;

import org.jlab.rec.dc.cross.Cross;

import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.Particle;
/*My track class simply extends the R3Cross with momentum-charge information.
 * 
 * Basically, this object is created for a Track with a R3 Cross, and stores the R3 cross info and momentum-charge info (from the track)
 * 
/* Also add info regarding the fact this track is associated with a generated particle or not
 * (If not, this is a "fake" track)
 */
public class TriggerTrack extends Cross{


    /**
     * 
     */
    
    private static final long serialVersionUID = 1L;
    private byte charge;
    private Vector3 momentum;
    private Particle genParticle=new Particle();
    private Boolean isAssociatedToGenParticle=false;
    
    public byte getCharge() {
        return charge;
    }

    public void setCharge(byte charge) {
        this.charge = charge;
    }

    public Vector3 getMomentum() {
        return momentum;
    }

    public void setMomentum(Vector3 momentum) {
        this.momentum = momentum;
    }
    
    public TriggerTrack(int sector) {
        super(sector,3,0);
    }
    
    public Particle getGenParticle(){
        return genParticle;
    }
    
    public void setGenParticle(Particle particle){
        genParticle.copy(particle);
        this.isAssociatedToGenParticle=true;
    }
    
    public Boolean isAssociatedToGenParticle(){
        return isAssociatedToGenParticle;
    }
    
    
}
