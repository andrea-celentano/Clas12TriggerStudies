package clas12;

import org.jlab.geom.prim.Point3D;

/*Since the AHit class in TOF package is quite complicate (and I do not need all the features), and the Hit class inherits from hit, it is better to create my own simple c
 * class
 */
public class RawTOFHit  implements HitWithDCPositionEnergyTimeInfo {
    private int _Sector;
    private int _Panel;
    private int _Paddle;
    private int _Id;
    private float _EnergyL;
    private float _EnergyR;
    private float _TimeL;
    private float _TimeR;
    private float _TimeAVG;
    
    private boolean isMatchedToR3Segments;
    private boolean isMatchedToR3CrossProjection;
    private double distanceR3CrossProjection;

    public RawTOFHit(int _Sector, int _Panel, int _Paddle, int _Id, float _EnergyL, float _EnergyR, float _TimeL, float _TimeR) {
        super();
        this._Panel = _Panel;
        this._Sector = _Sector;
        this._Paddle = _Paddle;
        this._EnergyL = _EnergyL;
        this._EnergyR = _EnergyR;
        this._TimeL = _TimeL;
        this._TimeR = _TimeR;
        this._Id = _Id;
        this._TimeAVG = (this._TimeL + this._TimeR) / 2;
    }

    public int get_Panel() {
        return _Panel;
    }

    public void set_Panel(int _Panel) {
        this._Panel = _Panel;
    }

    public int get_Sector() {
        return _Sector;
    }

    public void set_Sector(int _Sector) {
        this._Sector = _Sector;
    }

    public int get_Paddle() {
        return _Paddle;
    }

    public void set_Paddle(int _Paddle) {
        this._Paddle = _Paddle;
    }

    public int get_Id() {
        return _Id;
    }

    public void set_Id(int _Id) {
        this._Id = _Id;
    }

    public float get_EnergyL() {
        return _EnergyL;
    }

    public void set_EnergyL(float _EnergyL) {
        this._EnergyL = _EnergyL;
    }

    public float get_EnergyR() {
        return _EnergyR;
    }

    public void set_EnergyR(float _EnergyR) {
        this._EnergyR = _EnergyR;
    }

    public float get_TimeL() {
        return _TimeL;
    }

    public void set_TimeL(float _TimeL) {
        this._TimeL = _TimeL;
    }

    public float get_TimeR() {
        return _TimeR;
    }

    public void set_TimeR(float _TimeR) {
        this._TimeR = _TimeR;
    }

    public float get_TimeAVG() {
        return _TimeAVG;
    }

    public void set_TimeAVG(float _TimeAVG) {
        this._TimeAVG = _TimeAVG;
    }
    
    
    @Override
    public boolean isMatchedToR3Segments() {
        return isMatchedToR3Segments;
    }

    @Override
    public boolean isMatchedToR3CrossProjection() {
        return isMatchedToR3CrossProjection;
    }

    @Override
    public double distanceR3CrossProjection() {
        return distanceR3CrossProjection;
    }

    @Override
    public void setIsMatchedToR3Segments(boolean val) {
        if (!val) {
            isMatchedToR3CrossProjection = false;
            isMatchedToR3Segments = false;
        } else {
            isMatchedToR3CrossProjection = true;

        }
    }

    @Override
    public void setIsMatchedToR3CrossProjection(boolean val) {
        if (val) {
            isMatchedToR3CrossProjection = true;
            isMatchedToR3Segments = true;
        } else {
            isMatchedToR3CrossProjection = false;
        }
    }

    @Override
    public void setDistanceR3CrossProjection(double val) {
        this.distanceR3CrossProjection=val;

    }

    @Override
    public Point3D getPosition() {
        return null;
    }

    @Override
    public double getTime() {   
        return 0;
    }

    @Override
    public double getEnergy() {     
        return 0;
    }
    
    
    
}
