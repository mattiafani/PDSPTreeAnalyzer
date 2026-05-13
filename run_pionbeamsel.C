#include <map>

#include "./Analyzers/pionbeamsel.h"

void run_pionbeamsel() {
    // ─── USER CONFIGURATION ──────────────────────────────────────────
    const TString runMode = "MC";  // "MC"
    // const TString runMode = "Data";   // "Data"
    const double beamMomentum = 0.5;  // GeV/c: 0.5, 1.0, or 2.0
    const int maxEvents = -1;         // -1 = all
    // ─────────────────────────────────────────────────────────────────

    std::map<TString, std::map<double, TString>> fileMap = {
        {"MC", {{0.5, "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/22/59/77/18/PDSPProd4a_MC_0.5GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_04.root"}, {1.0, "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/18/80/01/67/PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"}, {2.0, "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/20/91/32/85/PDSPProd4a_MC_2GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"}}},
        {"Data", {{0.5, "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/58/25/PDSPProd4_data_0.5GeV_reco2_ntuple_v09_41_00_04.root"}, {1.0, "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/52/19/PDSPProd4_data_1GeV_reco2_ntuple_v09_41_00_04.root"}, {2.0, "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/54/29/PDSPProd4_data_2GeV_reco2_ntuple_v09_42_03_01.root"}}}};

    if (fileMap.find(runMode) == fileMap.end() ||
        fileMap[runMode].find(beamMomentum) == fileMap[runMode].end()) {
        std::cerr << "ERROR: no file registered for runMode=" << runMode
                  << " beamMomentum=" << beamMomentum << " GeV/c" << std::endl;
        return;
    }

    gSystem->Load("./lib/libDataFormats.so");
    gSystem->Load("./lib/libAnalyzerTools.so");
    gSystem->Load("./lib/libAnalyzers.so");

    pionbeamsel m;
    m.MaxEvent = maxEvents;
    m.LogEvery = 1000;
    m.MCSample = runMode;
    m.Beam_Momentum = beamMomentum;
    m.SetTreeName();
    m.AddFile(fileMap[runMode][beamMomentum]);

    TString outfile = TString::Format("hists_%s_%.1fGeV_pionbeamsel.root",
                                      runMode.Data(), beamMomentum);
    m.SetOutfilePath(outfile);
    m.Init();
    m.initializeAnalyzer();
    m.initializeAnalyzerTools();
    m.SwitchToTempDir();
    m.Loop();
    m.WriteHist();
}
