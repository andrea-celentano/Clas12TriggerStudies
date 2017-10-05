package clas12;

import java.util.ArrayList;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;

public class GuiClass {

    double deltaT; // how many ns per event?

    /* Here goes constants */
    static final int nCountersFTOF1A = 23;
    static final int nCountersFTOF1B = 62;
    static final int nCountersFTOF2 = 5;

    /* Here goes histograms */
    ArrayList<H1F> h1_rateFTOF1A;
    ArrayList<H1F> h1_rateFTOF1B;
    ArrayList<H1F> h1_rateFTOF2;

    ArrayList<H1F> h1_minPCAL_Distance;
    ArrayList<H1F> h1_allPCAL_E;
    ArrayList<H1F> h1_matchedPCAL_E;

    ArrayList<H1F> h1_minTOF1A_Distance;
    ArrayList<H1F> h1_allTOF1A_E;
    ArrayList<H1F> h1_matchedTOF1A_E;

    ArrayList<H1F> h1_minTOF1B_Distance;
    ArrayList<H1F> h1_allTOF1B_E;
    ArrayList<H1F> h1_matchedTOF1B_E;

    ArrayList<H1F> h1_minTOF2_Distance;
    ArrayList<H1F> h1_allTOF2_E;
    ArrayList<H1F> h1_matchedTOF2_E;

    H2F h2_ThetaPhiAllGEN;
    H2F h2_ThetaPhiAllREC;
    H2F h2_ThetaPhiAllTRIGGER;
    H2F h2_ThetaPhiAllTRIGGER2;
    H2F h2_ThetaPhiAllEFF;
    H2F h2_ThetaPhiAllEFF2;

    H2F h2_ThetaPhiQPGEN;
    H2F h2_ThetaPhiQPREC;
    H2F h2_ThetaPhiQPTRIGGER;
    H2F h2_ThetaPhiQPTRIGGER2;
    H2F h2_ThetaPhiQPEFF;
    H2F h2_ThetaPhiQPEFF2;

    H2F h2_ThetaPhiQMGEN;
    H2F h2_ThetaPhiQMREC;
    H2F h2_ThetaPhiQMTRIGGER;
    H2F h2_ThetaPhiQMTRIGGER2;
    H2F h2_ThetaPhiQMEFF;
    H2F h2_ThetaPhiQMEFF2;

    ArrayList<H2F> h2_ThetaPhiAllGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiAllREC_vsP;
    ArrayList<H2F> h2_ThetaPhiAllTRIGGER_vsP;
    ArrayList<H2F> h2_ThetaPhiAllTRIGGER2_vsP;
    ArrayList<H2F> h2_ThetaPhiAllEFF_vsP;
    ArrayList<H2F> h2_ThetaPhiAllEFF2_vsP;

    ArrayList<H2F> h2_ThetaPhiQPGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiQPREC_vsP;
    ArrayList<H2F> h2_ThetaPhiQPTRIGGER_vsP;
    ArrayList<H2F> h2_ThetaPhiQPTRIGGER2_vsP;
    ArrayList<H2F> h2_ThetaPhiQPEFF_vsP;
    ArrayList<H2F> h2_ThetaPhiQPEFF2_vsP;

    ArrayList<H2F> h2_ThetaPhiQMGEN_vsP;
    ArrayList<H2F> h2_ThetaPhiQMREC_vsP;
    ArrayList<H2F> h2_ThetaPhiQMTRIGGER_vsP;
    ArrayList<H2F> h2_ThetaPhiQMTRIGGER2_vsP;
    ArrayList<H2F> h2_ThetaPhiQMEFF_vsP;
    ArrayList<H2F> h2_ThetaPhiQMEFF2_vsP;

    H1F h1_vsDistanceAllREC;
    H1F h1_vsDistanceAllTRIGGER;
    H1F h1_vsDistanceAllTRIGGER2;
    H1F h1_vsDistanceAllEFF;
    H1F h1_vsDistanceAllEFF2;

    H1F h1_vsDistanceQPREC;
    H1F h1_vsDistanceQPTRIGGER;
    H1F h1_vsDistanceQPTRIGGER2;
    H1F h1_vsDistanceQPEFF;
    H1F h1_vsDistanceQPEFF2;

    H1F h1_vsDistanceQMREC;
    H1F h1_vsDistanceQMTRIGGER;
    H1F h1_vsDistanceQMTRIGGER2;
    H1F h1_vsDistanceQMEFF;
    H1F h1_vsDistanceQMEFF2;

    ArrayList<H1F> h1_vsDistanceAllREC_vsP;
    ArrayList<H1F> h1_vsDistanceAllTRIGGER_vsP;
    ArrayList<H1F> h1_vsDistanceAllTRIGGER2_vsP;
    ArrayList<H1F> h1_vsDistanceAllEFF_vsP;
    ArrayList<H1F> h1_vsDistanceAllEFF2_vsP;

    ArrayList<H1F> h1_vsDistanceQPREC_vsP;
    ArrayList<H1F> h1_vsDistanceQPTRIGGER_vsP;
    ArrayList<H1F> h1_vsDistanceQPTRIGGER2_vsP;
    ArrayList<H1F> h1_vsDistanceQPEFF_vsP;
    ArrayList<H1F> h1_vsDistanceQPEFF2_vsP;

    ArrayList<H1F> h1_vsDistanceQMREC_vsP;
    ArrayList<H1F> h1_vsDistanceQMTRIGGER_vsP;
    ArrayList<H1F> h1_vsDistanceQMTRIGGER2_vsP;
    ArrayList<H1F> h1_vsDistanceQMEFF_vsP;
    ArrayList<H1F> h1_vsDistanceQMEFF2_vsP;

    H1F h1_vsDistance1TRIGGER;
    H1F h1_vsDistance1TRIGGER2;
    H1F h1_vsDistance1EFF;
    H1F h1_vsDistance1EFF2;

    H1F h1_vsDistance2TRIGGER;
    H1F h1_vsDistance2TRIGGER2;
    H1F h1_vsDistance2EFF;
    H1F h1_vsDistance2EFF2;

    H1F h1_vsDistance3TRIGGER;
    H1F h1_vsDistance3TRIGGER2;
    H1F h1_vsDistance3EFF;
    H1F h1_vsDistance3EFF2;

    H2F h2_FTOF2EnergyAll_LR;
    H2F h2_FTOF2EnergyMatched_LR;

    H2F h2_FTOF1BEnergyAll_LR;
    H2F h2_FTOF1BEnergyMatched_LR;

    H2F h2_FTOF1AEnergyAll_LR;
    H2F h2_FTOF1AEnergyMatched_LR;

    H2F h2_FTOF1A_1B_match;

    ArrayList<H1F> allH1F;
    ArrayList<H2F> allH2F;

    ArrayList<Double> Parray;

    AnalysisClass analysisClass;

    public double getDeltaT() {
        return deltaT;
    }

    public void setDeltaT(double deltaT) {
        this.deltaT = deltaT;
    }

    public GuiClass(AnalysisClass ana) {
        this.analysisClass = ana;
    }

    public ArrayList<Double> getParray() {
        return Parray;
    }

    public void setParray(ArrayList<Double> parray) {
        Parray = parray;
    }

    public void setupHistograms() {

        this.allH1F = new ArrayList<H1F>();
        this.allH2F = new ArrayList<H2F>();

        h1_rateFTOF1A = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_rateFTOF1A.add(new H1F("h1_rateFTOF1A_" + isector, "h1_rateFTOF1A_" + isector, GuiClass.nCountersFTOF1A, 0.5, GuiClass.nCountersFTOF1A + .5));
            allH1F.add(h1_rateFTOF1A.get(isector - 1));
        }
        h1_rateFTOF1B = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_rateFTOF1B.add(new H1F("h1_rateFTOF1B_" + isector, "h1_rateFTOF1B_" + isector, GuiClass.nCountersFTOF1B, 0.5, GuiClass.nCountersFTOF1B + .5));
            allH1F.add(h1_rateFTOF1B.get(isector - 1));
        }
        h1_rateFTOF2 = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_rateFTOF2.add(new H1F("h1_rateFTOF2_" + isector, "h1_rateFTOF2_" + isector, GuiClass.nCountersFTOF2, 0.5, GuiClass.nCountersFTOF2 + .5));
            allH1F.add(h1_rateFTOF2.get(isector - 1));
        }

        h1_minPCAL_Distance = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_minPCAL_Distance.add(new H1F("h1_minPCAL_Distance_" + isector, "h1_minPCAL_Distance_" + isector + ";distance - cm", 100, 0, 500));
            allH1F.add(h1_minPCAL_Distance.get(isector - 1));
        }

        h1_allPCAL_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allPCAL_E.add(new H1F("h1_allPCAL_E_" + isector, "h1_allPCAL_E_" + isector + ";Energy (GeV)", 100, 0, .3));
            allH1F.add(h1_allPCAL_E.get(isector - 1));
        }

        h1_matchedPCAL_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_matchedPCAL_E.add(new H1F("h1_matchedPCAL_E_" + isector, "h1_matchedPCAL_E_" + isector + ";Energy (GeV)", 100, 0, .3));
            allH1F.add(h1_matchedPCAL_E.get(isector - 1));
        }

        h1_minTOF1A_Distance = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_minTOF1A_Distance.add(new H1F("h1_minTOF1A_Distance_" + isector, "h1_minTOF1A_Distance_" + isector + ";distance - cm", 100, 0, 500));
            allH1F.add(h1_minTOF1A_Distance.get(isector - 1));
        }

        h1_allTOF1A_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allTOF1A_E.add(new H1F("h1_allTOF1A_E_" + isector, "h1_allTOF1A_E_" + isector + ";Energy (GeV)", 100, 0, 30));
            allH1F.add(h1_allTOF1A_E.get(isector - 1));
        }

        h1_matchedTOF1A_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_matchedTOF1A_E.add(new H1F("h1_matchedTOF1A_E_" + isector, "h1_matchedTOF1A_E_" + isector + ";Energy (GeV)", 100, 0, 30));
            allH1F.add(h1_matchedTOF1A_E.get(isector - 1));
        }

        h1_minTOF1B_Distance = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_minTOF1B_Distance.add(new H1F("h1_minTOF1B_Distance_" + isector, "h1_minTOF1B_Distance_" + isector + ";distance - cm", 100, 0, 500));
            allH1F.add(h1_minTOF1B_Distance.get(isector - 1));
        }

        h1_allTOF1B_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allTOF1B_E.add(new H1F("h1_allTOF1B_E_" + isector, "h1_allTOF1B_E_" + isector + ";Energy (GeV)", 100, 0, 30));
            allH1F.add(h1_allTOF1B_E.get(isector - 1));
        }

        h1_matchedTOF1B_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_matchedTOF1B_E.add(new H1F("h1_matchedTOF1B_E_" + isector, "h1_matchedTOF1B_E_" + isector + ";Energy (GeV)", 100, 0, 30));
            allH1F.add(h1_matchedTOF1B_E.get(isector - 1));
        }

        h1_minTOF2_Distance = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_minTOF2_Distance.add(new H1F("h1_minTOF2_Distance_" + isector, "h1_minTOF2_Distance_" + isector + ";distance - cm", 100, 0, 500));
            allH1F.add(h1_minTOF2_Distance.get(isector - 1));
        }

        h1_allTOF2_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= 6; isector++) {
            h1_allTOF2_E.add(new H1F("h1_allTOF2_E_" + isector, "h1_allTOF2_E_" + isector + ";Energy (GeV)", 100, 0, 30));
            allH1F.add(h1_allTOF2_E.get(isector - 1));
        }

        h1_matchedTOF2_E = new ArrayList<H1F>();
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_matchedTOF2_E.add(new H1F("h1_matchedTOF2_E_" + isector, "h1_matchedTOF2_E_" + isector + ";Energy (GeV)", 100, 0, 30));
            allH1F.add(h1_matchedTOF2_E.get(isector - 1));
        }

        h2_ThetaPhiAllGEN = new H2F("h2_ThetaPhiAllGEN", "h2_ThetaPhiAllGEN", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllGEN);
        h2_ThetaPhiAllREC = new H2F("h2_ThetaPhiAllREC", "h2_ThetaPhiAllREC", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllREC);
        h2_ThetaPhiAllTRIGGER = new H2F("h2_ThetaPhiAllTRIGGER", "h2_ThetaPhiAllTRIGGER", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllTRIGGER);
        h2_ThetaPhiAllTRIGGER2 = new H2F("h2_ThetaPhiAllTRIGGER2", "h2_ThetaPhiAllTRIGGER2", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiAllTRIGGER2);

        h2_ThetaPhiQPGEN = new H2F("h2_ThetaPhiQPGEN", "h2_ThetaPhiQPGEN", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPGEN);
        h2_ThetaPhiQPREC = new H2F("h2_ThetaPhiQPREC", "h2_ThetaPhiQPREC", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPREC);
        h2_ThetaPhiQPTRIGGER = new H2F("h2_ThetaPhiQPTRIGGER", "h2_ThetaPhiQPTRIGGER", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPTRIGGER);
        h2_ThetaPhiQPTRIGGER2 = new H2F("h2_ThetaPhiQPTRIGGER2", "h2_ThetaPhiQPTRIGGER2", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQPTRIGGER2);

        h2_ThetaPhiQMGEN = new H2F("h2_ThetaPhiQMGEN", "h2_ThetaPhiQMGEN", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMGEN);
        h2_ThetaPhiQMREC = new H2F("h2_ThetaPhiQMREC", "h2_ThetaPhiQMREC", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMREC);
        h2_ThetaPhiQMTRIGGER = new H2F("h2_ThetaPhiQMTRIGGER", "h2_ThetaPhiQMTRIGGER", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMTRIGGER);
        h2_ThetaPhiQMTRIGGER2 = new H2F("h2_ThetaPhiQMTRIGGER2", "h2_ThetaPhiQMTRIGGER2", 100, -180., 180., 100, 0, 60.);
        allH2F.add(h2_ThetaPhiQMTRIGGER2);

        h1_vsDistanceAllREC = new H1F("h1_vsDistanceAllREC", "h1_vsDistanceAllREC", 400, 0., 400.);
        allH1F.add(h1_vsDistanceAllREC);
        h1_vsDistanceQPREC = new H1F("h1_vsDistanceQPREC", "h1_vsDistanceQPREC", 400, 0., 400.);
        allH1F.add(h1_vsDistanceQPREC);
        h1_vsDistanceQMREC = new H1F("h1_vsDistanceQMREC", "h1_vsDistanceQMREC", 400, 0., 400.);
        allH1F.add(h1_vsDistanceQMREC);

        h1_vsDistance1TRIGGER = new H1F("h1_vsDistance1TRIGGER", "h1_vsDistance3TRIGGER", 400, 0., 400.);
        allH1F.add(h1_vsDistance1TRIGGER);
        h1_vsDistance2TRIGGER = new H1F("h1_vsDistance2TRIGGER", "h1_vsDistance2TRIGGER", 400, 0., 400.);
        allH1F.add(h1_vsDistance2TRIGGER);
        h1_vsDistance3TRIGGER = new H1F("h1_vsDistance3TRIGGER", "h1_vsDistance3TRIGGER", 400, 0., 400.);
        allH1F.add(h1_vsDistance3TRIGGER);

        h1_vsDistance1TRIGGER2 = new H1F("h1_vsDistance1TRIGGER2", "h1_vsDistance3TRIGGER2", 400, 0., 400.);
        allH1F.add(h1_vsDistance1TRIGGER2);
        h1_vsDistance2TRIGGER2 = new H1F("h1_vsDistance2TRIGGER2", "h1_vsDistance2TRIGGER2", 400, 0., 400.);
        allH1F.add(h1_vsDistance2TRIGGER2);
        h1_vsDistance3TRIGGER2 = new H1F("h1_vsDistance3TRIGGER2", "h1_vsDistance3TRIGGER2", 400, 0., 400.);
        allH1F.add(h1_vsDistance3TRIGGER2);

        h1_vsDistanceAllREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllTRIGGER_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllTRIGGER2_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllEFF_vsP = new ArrayList<H1F>();
        h1_vsDistanceAllEFF2_vsP = new ArrayList<H1F>();

        h1_vsDistanceQPREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPTRIGGER_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPTRIGGER2_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPEFF_vsP = new ArrayList<H1F>();
        h1_vsDistanceQPEFF2_vsP = new ArrayList<H1F>();

        h1_vsDistanceQMREC_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMTRIGGER_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMTRIGGER2_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMEFF_vsP = new ArrayList<H1F>();
        h1_vsDistanceQMEFF2_vsP = new ArrayList<H1F>();

        h2_ThetaPhiAllGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllTRIGGER_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllTRIGGER2_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllEFF_vsP = new ArrayList<H2F>();
        h2_ThetaPhiAllEFF2_vsP = new ArrayList<H2F>();

        h2_ThetaPhiQPGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPTRIGGER_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPTRIGGER2_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPEFF_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQPEFF2_vsP = new ArrayList<H2F>();

        h2_ThetaPhiQMGEN_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMREC_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMTRIGGER_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMTRIGGER2_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMEFF_vsP = new ArrayList<H2F>();
        h2_ThetaPhiQMEFF2_vsP = new ArrayList<H2F>();

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {

            h1_vsDistanceAllREC_vsP.add(new H1F("h1_vsDistanceAllREC_" + ibin, "h1_vsDistanceAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPREC_vsP.add(new H1F("h1_vsDistanceQPREC_" + ibin, "h1_vsDistanceQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMREC_vsP.add(new H1F("h1_vsDistanceQMREC_" + ibin, "h1_vsDistanceQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            h1_vsDistanceAllTRIGGER_vsP.add(new H1F("h1_vsDistanceAllTRIGGER_" + ibin, "h1_vsDistanceAllTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPTRIGGER_vsP.add(new H1F("h1_vsDistanceQPTRIGGER_" + ibin, "h1_vsDistanceQPTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMTRIGGER_vsP.add(new H1F("h1_vsDistanceQMTRIGGER_" + ibin, "h1_vsDistanceQMTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            h1_vsDistanceAllTRIGGER2_vsP.add(new H1F("h1_vsDistanceAllTRIGGER2_" + ibin, "h1_vsDistanceAllTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQPTRIGGER2_vsP.add(new H1F("h1_vsDistanceQPTRIGGER2_" + ibin, "h1_vsDistanceQPTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));
            h1_vsDistanceQMTRIGGER2_vsP.add(new H1F("h1_vsDistanceQMTRIGGER2_" + ibin, "h1_vsDistanceQMTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 400, 0., 400.));

            allH1F.add(h1_vsDistanceAllREC_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQPREC_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQMREC_vsP.get(ibin));
            allH1F.add(h1_vsDistanceAllTRIGGER_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQPTRIGGER_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQMTRIGGER_vsP.get(ibin));
            allH1F.add(h1_vsDistanceAllTRIGGER2_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQPTRIGGER2_vsP.get(ibin));
            allH1F.add(h1_vsDistanceQMTRIGGER2_vsP.get(ibin));

            h2_ThetaPhiAllGEN_vsP.add(new H2F("h2_ThetaPhiAllGEN_" + ibin, "h2_ThetaPhiAllGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiAllREC_vsP.add(new H2F("h2_ThetaPhiAllREC_" + ibin, "h2_ThetaPhiAllREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiAllTRIGGER_vsP.add(new H2F("h2_ThetaPhiAllTRIGGER_" + ibin, "h2_ThetaPhiAllTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiAllTRIGGER2_vsP.add(new H2F("h2_ThetaPhiAllTRIGGER2_" + ibin, "h2_ThetaPhiAllTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));

            allH2F.add(h2_ThetaPhiAllGEN_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiAllREC_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiAllTRIGGER_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiAllTRIGGER2_vsP.get(ibin));

            h2_ThetaPhiQPGEN_vsP.add(new H2F("h2_ThetaPhiQPGEN_" + ibin, "h2_ThetaPhiQPGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQPREC_vsP.add(new H2F("h2_ThetaPhiQPREC_" + ibin, "h2_ThetaPhiQPREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQPTRIGGER_vsP.add(new H2F("h2_ThetaPhiQPTRIGGER_" + ibin, "h2_ThetaPhiQPTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQPTRIGGER2_vsP.add(new H2F("h2_ThetaPhiQPTRIGGER2_" + ibin, "h2_ThetaPhiQPTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));

            allH2F.add(h2_ThetaPhiQPGEN_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQPREC_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQPTRIGGER_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQPTRIGGER2_vsP.get(ibin));

            h2_ThetaPhiQMGEN_vsP.add(new H2F("h2_ThetaPhiQMGEN_" + ibin, "h2_ThetaPhiQMGEN_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQMREC_vsP.add(new H2F("h2_ThetaPhiQMREC_" + ibin, "h2_ThetaPhiQMREC_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQMTRIGGER_vsP.add(new H2F("h2_ThetaPhiQMTRIGGER_" + ibin, "h2_ThetaPhiQMTRIGGER_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));
            h2_ThetaPhiQMTRIGGER2_vsP.add(new H2F("h2_ThetaPhiQMTRIGGER2_" + ibin, "h2_ThetaPhiQMTRIGGER2_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1), 100, -180., 180., 100, 0, 60.));

            allH2F.add(h2_ThetaPhiQMGEN_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQMREC_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQMTRIGGER_vsP.get(ibin));
            allH2F.add(h2_ThetaPhiQMTRIGGER2_vsP.get(ibin));

        }

        /* MIP: ~ 2 MeV / cm, TOF2 thick=5 cm -> MPI=10 MeV */
        h2_FTOF2EnergyAll_LR = new H2F("h2_FTOF2EnergyAll_LR", "h2_FTOF2EnergyAll_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF2EnergyAll_LR);

        h2_FTOF2EnergyMatched_LR = new H2F("h2_FTOF2EnergyMatched_LR", "h2_FTOF2EnergyMatched_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF2EnergyMatched_LR);

        h2_FTOF1BEnergyAll_LR = new H2F("h2_FTOF1BEnergyAll_LR", "h2_FTOF1BEnergyAll_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1BEnergyAll_LR);

        h2_FTOF1BEnergyMatched_LR = new H2F("h2_FTOF1BEnergyMatched_LR", "h2_FTOF1BEnergyMatched_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1BEnergyMatched_LR);

        h2_FTOF1AEnergyAll_LR = new H2F("h2_FTOF1AEnergyAll_LR", "h2_FTOF1AEnergyAll_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1AEnergyAll_LR);

        h2_FTOF1AEnergyMatched_LR = new H2F("h2_FTOF1AEnergyMatched_LR", "h2_FTOF1AEnergyMatched_LR", 100, 0, 30, 100, 0, 30);
        allH2F.add(h2_FTOF1AEnergyMatched_LR);

        h2_FTOF1A_1B_match = new H2F("h2_FTOF1A_1B_match", "h2_FTOF1A_1B_match", 25, -0.5, 24.5, 65, -0.5, 64.5);
        allH2F.add(h2_FTOF1A_1B_match);
    }

    public void showHistograms() {

        TCanvas trk1 = new TCanvas("trk1", 1600, 1000);
        trk1.divide(4, 3);
        trk1.cd(0);
        for (int isector = 1; isector <= 6; isector++) {
            h1_minPCAL_Distance.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_minPCAL_Distance.get(isector - 1));
            } else {
                trk1.draw(h1_minPCAL_Distance.get(isector - 1), "same");
            }
        }
        trk1.cd(4);
        for (int isector = 1; isector <= 6; isector++) {
            h1_allPCAL_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_allPCAL_E.get(isector - 1));
            } else {
                trk1.draw(h1_allPCAL_E.get(isector - 1), "same");
            }
        }
        trk1.cd(8);
        for (int isector = 1; isector <= 6; isector++) {
            h1_matchedPCAL_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_matchedPCAL_E.get(isector - 1));
            } else {
                trk1.draw(h1_matchedPCAL_E.get(isector - 1), "same");
            }
        }

        trk1.cd(1);
        for (int isector = 1; isector <= 6; isector++) {
            h1_minTOF1A_Distance.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_minTOF1A_Distance.get(isector - 1));
            } else {
                trk1.draw(h1_minTOF1A_Distance.get(isector - 1), "same");
            }
        }
        trk1.cd(5);
        for (int isector = 1; isector <= 6; isector++) {
            h1_allTOF1A_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_allTOF1A_E.get(isector - 1));
            } else {
                trk1.draw(h1_allTOF1A_E.get(isector - 1), "same");
            }
        }
        trk1.cd(9);
        for (int isector = 1; isector <= 6; isector++) {
            h1_matchedTOF1A_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_matchedTOF1A_E.get(isector - 1));
            } else {
                trk1.draw(h1_matchedTOF1A_E.get(isector - 1), "same");
            }
        }

        trk1.cd(2);
        for (int isector = 1; isector <= 6; isector++) {
            h1_minTOF1B_Distance.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_minTOF1B_Distance.get(isector - 1));
            } else {
                trk1.draw(h1_minTOF1B_Distance.get(isector - 1), "same");
            }
        }
        trk1.cd(6);
        for (int isector = 1; isector <= 6; isector++) {
            h1_allTOF1B_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_allTOF1B_E.get(isector - 1));
            } else {
                trk1.draw(h1_allTOF1B_E.get(isector - 1), "same");
            }
        }
        trk1.cd(10);
        for (int isector = 1; isector <= 6; isector++) {
            h1_matchedTOF1B_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_matchedTOF1B_E.get(isector - 1));
            } else {
                trk1.draw(h1_matchedTOF1B_E.get(isector - 1), "same");
            }
        }

        trk1.cd(3);
        for (int isector = 1; isector <= 6; isector++) {
            h1_minTOF2_Distance.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_minTOF2_Distance.get(isector - 1));
            } else {
                trk1.draw(h1_minTOF2_Distance.get(isector - 1), "same");
            }
        }
        trk1.cd(7);
        for (int isector = 1; isector <= 6; isector++) {
            h1_allTOF2_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_allTOF2_E.get(isector - 1));
            } else {
                trk1.draw(h1_allTOF2_E.get(isector - 1), "same");
            }
        }
        trk1.cd(11);
        for (int isector = 1; isector <= 6; isector++) {
            h1_matchedTOF2_E.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                trk1.draw(h1_matchedTOF2_E.get(isector - 1));
            } else {
                trk1.draw(h1_matchedTOF2_E.get(isector - 1), "same");
            }
        }

        TCanvas trkAllP = new TCanvas("allP", 1600, 1000);
        trkAllP.divide(5, 3);
        trkAllP.cd(0);
        trkAllP.draw(h1_vsDistanceAllEFF);
        h1_vsDistanceQPEFF.setLineColor(2);
        trkAllP.draw(h1_vsDistanceQPEFF, "same");
        h1_vsDistanceQMEFF.setLineColor(3);
        trkAllP.draw(h1_vsDistanceQMEFF, "same");

        trkAllP.cd(5);
        trkAllP.draw(h2_ThetaPhiQPGEN, "colz");

        trkAllP.cd(6);
        trkAllP.draw(h2_ThetaPhiQPREC, "colz");

        trkAllP.cd(7);
        trkAllP.draw(h2_ThetaPhiQPTRIGGER, "colz");

        trkAllP.cd(8);
        trkAllP.draw(h2_ThetaPhiQPTRIGGER2, "colz");

        trkAllP.cd(9);
        trkAllP.getCanvas().getPad(9).getAxisZ().setLog(false);
        trkAllP.draw(h2_ThetaPhiQPEFF, "colz");

        trkAllP.cd(10);
        trkAllP.draw(h2_ThetaPhiQMGEN, "colz");

        trkAllP.cd(11);
        trkAllP.draw(h2_ThetaPhiQMREC, "colz");

        trkAllP.cd(12);
        trkAllP.draw(h2_ThetaPhiQMTRIGGER, "colz");

        trkAllP.cd(13);
        trkAllP.draw(h2_ThetaPhiQMTRIGGER2, "colz");

        trkAllP.cd(14);
        trkAllP.getCanvas().getPad(14).getAxisZ().setLog(false);
        trkAllP.draw(h2_ThetaPhiQMEFF, "colz");

        ArrayList<TCanvas> TCanvasArray = new ArrayList<TCanvas>();
        for (int ii = 0; ii < Parray.size() - 1; ii++) {
            TCanvasArray.add(new TCanvas("P_" + Parray.get(ii) + "_" + Parray.get(ii + 1), 1600, 1000));
            TCanvasArray.get(ii).divide(5, 3);

            TCanvasArray.get(ii).cd(0);
            TCanvasArray.get(ii).draw(h1_vsDistanceAllEFF_vsP.get(ii));
            h1_vsDistanceQPEFF_vsP.get(ii).setLineColor(2);
            TCanvasArray.get(ii).draw(h1_vsDistanceQPEFF_vsP.get(ii), "same");
            h1_vsDistanceQMEFF_vsP.get(ii).setLineColor(3);
            TCanvasArray.get(ii).draw(h1_vsDistanceQMEFF_vsP.get(ii), "same");

            TCanvasArray.get(ii).cd(5);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPGEN_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(6);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPREC_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(7);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPTRIGGER_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(8);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPTRIGGER2_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(9);
            TCanvasArray.get(ii).getCanvas().getPad(9).getAxisZ().setLog(false);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQPEFF_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(10);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMGEN_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(11);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMREC_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(12);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMTRIGGER_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(13);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMTRIGGER2_vsP.get(ii), "colz");

            TCanvasArray.get(ii).cd(14);
            TCanvasArray.get(ii).getCanvas().getPad(14).getAxisZ().setLog(false);
            TCanvasArray.get(ii).draw(h2_ThetaPhiQMEFF_vsP.get(ii), "colz");

        }

        TCanvas ceff = new TCanvas("ceff", 1600, 1000);
        ceff.divide(2, 3);

        ceff.cd(0);
        ceff.draw(h1_vsDistance1EFF);
        h1_vsDistance2EFF.setLineColor(2);
        ceff.draw(h1_vsDistance2EFF, "same");
        h1_vsDistance3EFF.setLineColor(3);
        ceff.draw(h1_vsDistance3EFF, "same");

        ceff.cd(1);
        ceff.draw(h1_vsDistance1TRIGGER);
        h1_vsDistance2TRIGGER.setLineColor(2);
        ceff.draw(h1_vsDistance2TRIGGER, "same");
        h1_vsDistance3TRIGGER.setLineColor(3);
        ceff.draw(h1_vsDistance3TRIGGER, "same");

        ceff.cd(2);
        ceff.draw(h1_vsDistance1EFF2);
        h1_vsDistance2EFF2.setLineColor(2);
        ceff.draw(h1_vsDistance2EFF2, "same");
        h1_vsDistance3EFF2.setLineColor(3);
        ceff.draw(h1_vsDistance3EFF2, "same");

        ceff.cd(3);
        ceff.draw(h1_vsDistance1TRIGGER2);
        h1_vsDistance2TRIGGER2.setLineColor(2);
        ceff.draw(h1_vsDistance2TRIGGER2, "same");
        h1_vsDistance3TRIGGER2.setLineColor(3);
        ceff.draw(h1_vsDistance3TRIGGER2, "same");

        TCanvas ceff2 = new TCanvas("ceff2", 1600, 1000);
        ceff2.divide(Parray.size() - 1, 1);
        for (int ii = 0; ii < Parray.size() - 1; ii++) {
            ceff2.cd(ii);
            ceff2.draw(h2_ThetaPhiAllEFF_vsP.get(ii), "colz");
        }

        TCanvas cTOF = new TCanvas("cTOF2", 1600, 1600);
        cTOF.divide(3, 3);

        cTOF.cd(0);
        cTOF.getCanvas().getPad(0).getAxisZ().setLog(false);
        cTOF.draw(h2_FTOF1AEnergyAll_LR, "colz");

        cTOF.cd(1);
        cTOF.getCanvas().getPad(1).getAxisZ().setLog(false);
        cTOF.draw(h2_FTOF1BEnergyAll_LR, "colz");

        cTOF.cd(2);
        cTOF.getCanvas().getPad(2).getAxisZ().setLog(false);
        cTOF.draw(h2_FTOF2EnergyAll_LR, "colz");

        cTOF.cd(3);
        cTOF.getCanvas().getPad(3).getAxisZ().setLog(false);
        cTOF.draw(h2_FTOF1AEnergyMatched_LR, "colz");

        cTOF.cd(4);
        cTOF.getCanvas().getPad(4).getAxisZ().setLog(false);
        cTOF.draw(h2_FTOF1BEnergyMatched_LR, "colz");

        cTOF.cd(5);
        cTOF.getCanvas().getPad(5).getAxisZ().setLog(false);
        cTOF.draw(h2_FTOF2EnergyMatched_LR, "colz");
        // cTOF2.draw(h2tmp,"colz");

        cTOF.cd(6);
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_rateFTOF1A.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                cTOF.draw(h1_rateFTOF1A.get(isector - 1));
            } else {
                cTOF.draw(h1_rateFTOF1A.get(isector - 1), "same");
            }
        }
        cTOF.cd(7);
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_rateFTOF1B.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                cTOF.draw(h1_rateFTOF1B.get(isector - 1));
            } else {
                cTOF.draw(h1_rateFTOF1B.get(isector - 1), "same");
            }
        }
        cTOF.cd(8);
        for (int isector = 1; isector <= AnalysisClass.nSectors_CLAS12; isector++) {
            h1_rateFTOF2.get(isector - 1).setLineColor(isector);
            if (isector == 1) {
                cTOF.draw(h1_rateFTOF2.get(isector - 1));
            } else {
                cTOF.draw(h1_rateFTOF2.get(isector - 1), "same");
            }
        }

    }

    public void processHistograms() {

        int nEvents = analysisClass.nevent;
        int nTracksWithR3Cross = analysisClass.nTracksMatchedToGen;
        int nTracksQPWithR3Cross = analysisClass.nTracksQPMatchedToGen;
        int nTracksQMWithR3Cross = analysisClass.nTracksQMMatchedToGen;
        double T = nEvents * this.deltaT; // in s
        T *= 1E3; // in ms, to have kHz rate

        System.out.println("Process Histogrqms: time (ms) is " + T);

        ArrayList<Integer> nTracks_vsP = analysisClass.nTracksMatchedToGen_vsP;
        ArrayList<Integer> nTracksQP_vsP = analysisClass.nTracksMatchedToGenQP_vsP;
        ArrayList<Integer> nTracksQM_vsP = analysisClass.nTracksMatchedToGenQM_vsP;

        h1_vsDistanceAllEFF = h1_vsDistanceAllREC.histClone("h1_vsDistanceAllEFF");
        h1_vsDistanceAllEFF.setTitle("h1_vsDistanceAllEFF");
        h1_vsDistanceAllEFF.divide(1. * nTracksWithR3Cross);

        h1_vsDistanceQPEFF = h1_vsDistanceQPREC.histClone("h1_vsDistanceQPEFF");
        h1_vsDistanceQPEFF.setTitle("h1_vsDistanceQPEFF");
        h1_vsDistanceQPEFF.divide(1. * nTracksQPWithR3Cross);

        h1_vsDistanceQMEFF = h1_vsDistanceQMREC.histClone("h1_vsDistanceQMEFF");
        h1_vsDistanceQMEFF.setTitle("h1_vsDistanceQMEFF");
        h1_vsDistanceQMEFF.divide(1. * nTracksQMWithR3Cross);

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {
            h1_vsDistanceAllEFF_vsP.add(h1_vsDistanceAllREC_vsP.get(ibin).histClone("h1_vsDistanceAllEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1)));
            h1_vsDistanceAllEFF_vsP.get(ibin).setTitle("h1_vsDistanceAllEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1));
            h1_vsDistanceAllEFF_vsP.get(ibin).divide(1. * nTracks_vsP.get(ibin));

            h1_vsDistanceQPEFF_vsP.add(h1_vsDistanceQPREC_vsP.get(ibin).histClone("h1_vsDistanceQPEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1)));
            h1_vsDistanceQPEFF_vsP.get(ibin).setTitle("h1_vsDistanceQPEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1));
            h1_vsDistanceQPEFF_vsP.get(ibin).divide(1. * nTracksQP_vsP.get(ibin));

            h1_vsDistanceQMEFF_vsP.add(h1_vsDistanceQMREC_vsP.get(ibin).histClone("h1_vsDistanceQMEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1)));
            h1_vsDistanceQMEFF_vsP.get(ibin).setTitle("h1_vsDistanceQMEFF_" + Parray.get(ibin) + "_" + Parray.get(ibin + 1));
            h1_vsDistanceQMEFF_vsP.get(ibin).divide(1. * nTracksQM_vsP.get(ibin));
        }

        h2_ThetaPhiAllEFF = H2F.divide(h2_ThetaPhiAllTRIGGER2, h2_ThetaPhiAllREC);
        h2_ThetaPhiQPEFF = H2F.divide(h2_ThetaPhiQPTRIGGER2, h2_ThetaPhiQPREC);
        h2_ThetaPhiQMEFF = H2F.divide(h2_ThetaPhiQMTRIGGER2, h2_ThetaPhiQMREC);

        for (int ibin = 0; ibin < Parray.size() - 1; ibin++) {
            h2_ThetaPhiAllEFF_vsP.add(H2F.divide(h2_ThetaPhiAllTRIGGER2_vsP.get(ibin), h2_ThetaPhiAllREC_vsP.get(ibin)));
            h2_ThetaPhiQPEFF_vsP.add(H2F.divide(h2_ThetaPhiQPTRIGGER2_vsP.get(ibin), h2_ThetaPhiQPREC_vsP.get(ibin)));
            h2_ThetaPhiQMEFF_vsP.add(H2F.divide(h2_ThetaPhiQMTRIGGER2_vsP.get(ibin), h2_ThetaPhiQMREC_vsP.get(ibin)));

            h2_ThetaPhiAllEFF_vsP.get(ibin).setTitle(Parray.get(ibin) + "-" + Parray.get(ibin + 1)+" GeV/c");
            h2_ThetaPhiQPEFF_vsP.get(ibin).setTitle(Parray.get(ibin) + "-" + Parray.get(ibin + 1)+" GeV/c");
            h2_ThetaPhiQMEFF_vsP.get(ibin).setTitle(Parray.get(ibin) + "-" + Parray.get(ibin + 1)+" GeV/c");

        }

        h1_vsDistance1EFF = h1_vsDistance1TRIGGER.histClone("h1_vsDistance1EFF");
        h1_vsDistance1EFF.setTitle("h1_vsDistance1EFF");
        h1_vsDistance1EFF.divide(1. * nEvents);

        h1_vsDistance2EFF = h1_vsDistance2TRIGGER.histClone("h1_vsDistance2EFF");
        h1_vsDistance2EFF.setTitle("h1_vsDistance2EFF");
        h1_vsDistance2EFF.divide(1. * nEvents);

        h1_vsDistance3EFF = h1_vsDistance3TRIGGER.histClone("h1_vsDistance3EFF");
        h1_vsDistance3EFF.setTitle("h1_vsDistance3EFF");
        h1_vsDistance3EFF.divide(1. * nEvents);

        h1_vsDistance1EFF2 = h1_vsDistance1TRIGGER2.histClone("h1_vsDistance1EFF2");
        h1_vsDistance1EFF2.setTitle("h1_vsDistance1EFF2");
        h1_vsDistance1EFF2.divide(1. * nEvents);

        h1_vsDistance2EFF2 = h1_vsDistance2TRIGGER2.histClone("h1_vsDistance2EFF2");
        h1_vsDistance2EFF2.setTitle("h1_vsDistance2EFF2");
        h1_vsDistance2EFF2.divide(1. * nEvents);

        h1_vsDistance3EFF2 = h1_vsDistance3TRIGGER2.histClone("h1_vsDistance3EFF2");
        h1_vsDistance3EFF2.setTitle("h1_vsDistance3EFF2");
        h1_vsDistance3EFF2.divide(1. * nEvents);

        for (int ii = 0; ii < AnalysisClass.nSectors_CLAS12; ii++) {
            h1_rateFTOF1A.get(ii).divide(T);
            h1_rateFTOF1B.get(ii).divide(T);
            h1_rateFTOF2.get(ii).divide(T);
        }

    }

    public H1F getHistogram1D(String name) {
        for (H1F h : this.allH1F) {
            if (name.equals(h.getName()) == true) {
                return h;
            }
        }
        System.out.println("histogram: " + name + " not found");
        return null;
    }

    public H2F getHistogram2D(String name) {
        for (H2F h : this.allH2F) {
            if (name.equals(h.getName()) == true) {
                return h;
            }
        }
        System.out.println("histogram: " + name + " not found");
        return null;
    }

}
