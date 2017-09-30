package clas12;
import org.jlab.rec.dc.cross.Cross;

/* Add info regarding the matching:
 * --> with calorimeter - in terms of the minimum distance to a good cluster (-1 if not matched)
 * --> with FTOF2 (y/n)
 * --> with FTOF1B (y/n)
*/

public class MatchedCross extends Cross {

    
  
    private static final long serialVersionUID = 1L;
    
    private Boolean isMatchedToFTOF1A=false;
    private Boolean isMatchedToFTOF1B=false;
    private Boolean isMatchedToFTOF2=false;
    private Boolean isMatchedToEC=false;

    
    private double minDistanceFTOF1A=-1;
    private double minDistanceFTOF1B=-1;
    private double minDistanceFTOF2=-1; 
    private double minDistanceEC=-1;
    
    public MatchedCross(int sector) {
        super(sector,3,0);
    }

    public Boolean isMatchedToFTOF1A() {
        return isMatchedToFTOF1A;
    }

    public void setIsMatchedToFTOF1A(Boolean isMatchedToFTOF1A) {
        this.isMatchedToFTOF1A = isMatchedToFTOF1A;
    }

    public Boolean isMatchedToFTOF1B() {
        return isMatchedToFTOF1B;
    }

    public void setIsMatchedToFTOF1B(Boolean isMatchedToFTOF1B) {
        this.isMatchedToFTOF1B = isMatchedToFTOF1B;
    }

    public Boolean isMatchedToFTOF2() {
        return isMatchedToFTOF2;
    }

    public void setIsMatchedToFTOF2(Boolean isMatchedToFTOF2) {
        this.isMatchedToFTOF2 = isMatchedToFTOF2;
    }

    public Boolean isMatchedToEC() {
        return isMatchedToEC;
    }

    public void setIsMatchedToEC(Boolean isMatchedToEC) {
        this.isMatchedToEC = isMatchedToEC;
    }

    public double getMinDistanceFTOF1A() {
        return minDistanceFTOF1A;
    }

    public void setMinDistanceFTOF1A(double minDistanceFTOF1A) {
        this.minDistanceFTOF1A = minDistanceFTOF1A;
    }

    public double getMinDistanceFTOF1B() {
        return minDistanceFTOF1B;
    }

    public void setMinDistanceFTOF1B(double minDistanceFTOF1B) {
        this.minDistanceFTOF1B = minDistanceFTOF1B;
    }

    public double getMinDistanceFTOF2() {
        return minDistanceFTOF2;
    }

    public void setMinDistanceFTOF2(double minDistanceFTOF2) {
        this.minDistanceFTOF2 = minDistanceFTOF2;
    }

    public double getMinDistanceEC() {
        return minDistanceEC;
    }

    public void setMinDistanceEC(double minDistanceEC) {
        this.minDistanceEC = minDistanceEC;
    }
    
   
}
