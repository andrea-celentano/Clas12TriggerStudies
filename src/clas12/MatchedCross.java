package clas12;
import org.jlab.rec.dc.cross.Cross;

/* Also add info regarding the matching:
 * --> with calorimeter - in terms of the minimum distance to a good cluster (-1 if not matched)
 * --> with FTOF2 (y/n)
*/

public class MatchedCross extends Cross {

    
  
    private static final long serialVersionUID = 1L;
    
    private Boolean isMatchedToFTOF2=false;
    private double minDistanceEC=-1;
    
    public MatchedCross(int sector) {
        super(sector,3,0);
    }
    
    public Boolean isMatchedToFTOF2() {
        return isMatchedToFTOF2;
    }

    public void setIsMatchedToFTOF2(Boolean isMatchedToFTOF2) {
        this.isMatchedToFTOF2 = isMatchedToFTOF2;
    }

    public double getMinDistanceEC() {
        return minDistanceEC;
    }

    public void setMinDistanceEC(double minDistanceEC) {
        this.minDistanceEC = minDistanceEC;
    }
}
