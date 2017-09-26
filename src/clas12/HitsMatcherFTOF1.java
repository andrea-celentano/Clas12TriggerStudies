package clas12;

import java.util.List;
/*This class handles the matching between FTOF1b and FTOF1a layers*/

public class HitsMatcherFTOF1 {

    private AnalysisClass analysisClass;
    private int sector;
    private double Tmin,Tcoinc;
    
    public HitsMatcherFTOF1(AnalysisClass ana_,int sect_,double Tcoinc_,double Tmin_) {
        this.analysisClass = ana_;
        this.sector = sect_;
        this.Tmin=Tmin_;
        this.Tcoinc=Tcoinc_;
    }
    
    /*This performs the matching.
     * Each list is made of hits having BOTH PMTs firing (I am using an energy thr over FTOF::rawhits in GEMC).
     * Simply return a boolean: has the sector a 1A/1B coincidence within the time window of length Tcoinc starting at Tmin? 
     */
    public boolean MatchHits(List<SimpleTOFHit> hitsFTOF1A,List<SimpleTOFHit> hitsFTOF1B){
    
        boolean ret=false;
        
        if (hitsFTOF1A.size()==0) return false;
        if (hitsFTOF1B.size()==0) return false;
        
        for (SimpleTOFHit hit1A : hitsFTOF1A){
            /*First, exclude hits not in the coinc. window*/            
            if (hit1A.get_TimeAVG()<Tmin) continue;
            if (hit1A.get_TimeAVG()>(Tmin+Tcoinc)) continue;
            
            for (SimpleTOFHit hit1B : hitsFTOF1B){
                /*Same for 1B, exclude hits not in the coinc. window*/            
                if (hit1B.get_TimeAVG()<Tmin) continue;
                if (hit1B.get_TimeAVG()>(Tmin+Tcoinc)) continue;
                
                if (this.MatchGeo(hit1A,hit1B)){
                    ret=true;
                    break;
                }
            }
        }
        
        return ret;
    }
    
    private boolean MatchGeo(SimpleTOFHit hit1A,SimpleTOFHit hit1B){
        boolean ret=true;
        return ret;
    }
    
}
