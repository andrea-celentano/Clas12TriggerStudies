package clas12;

public class SimpleDCSegment {
    private int sector;
    private int superlayer;
    private int layer;

    public SimpleDCSegment(int sector, int superlayer) {
        super();
        this.sector = sector;
        this.superlayer = superlayer;
        this.layer = (superlayer+1)/2;
    }

    public int getSector() {
        return sector;
    }

    public void setSector(int sector) {
        this.sector = sector;
    }

    public int getSuperlayer() {
        return superlayer;
    }

    public void setSuperlayer(int superlayer) {
        this.superlayer = superlayer;
    }

    public int getLayer() {
        return layer;
    }

    public void setLayer(int layer) {
        this.layer = layer;
    }
    
    
    
}
