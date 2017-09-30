package clas12;

public interface HitWithDCInfo {

    boolean isMatchedToR3Segments();
    boolean isMatchedToR3CrossProjection();
    double  distanceR3CrossProjection();
    
    void setIsMatchedToR3Segments(boolean val);
    void setIsMatchedToR3CrossProjection(boolean val);
    void setDistanceR3CrossProjection(double val);
}
