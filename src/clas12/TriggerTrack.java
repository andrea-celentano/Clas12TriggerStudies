package clas12;

import java.util.ArrayList;
import java.util.List;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.clas.physics.Vector3;

/*My track class simply extends the R3Cross with momentum-charge information*/
public class TriggerTrack extends Cross{


    private byte charge;
    private Vector3 momentum;


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
}
