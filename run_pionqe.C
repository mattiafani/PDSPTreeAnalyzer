#include "./Analyzers/pionqe0p5.h"

// ============================================================
// Run the pionqe0p5 analyzer over one sample (MC or Data).
// Called by run_pionqe(); not intended for direct use.
// ============================================================
void run_one_sample(TString sample) {
    if (sample != "MC" && sample != "Data") {
        printf("ERROR: sample must be \"MC\" or \"Data\" (got \"%s\")\n",
               sample.Data());
        return;
    }
    printf("\n========================================================\n");
    printf("  Running pionqe0p5 on %s sample\n", sample.Data());
    printf("========================================================\n");

    pionqe0p5 m;
    m.MaxEvent = -1;
    // m.MaxEvent = 3000;
    // m.NSkipEvent = 0;
    // m.MaxEvent = m.NSkipEvent + 1000;
    m.LogEvery = 1000;

    m.MCSample = sample;

    m.Beam_Momentum = 0.5;
    m.KE_ff_subt = 30.;
    m.SetTreeName();

    // ---------------- INPUT FILES ----------------
    // Historical / alternative samples (kept for reference):
    // m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/20/91/32/85/PDSPProd4a_MC_2GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"); // 2 GeV MC
    // m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/18/80/01/67/PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"); // 1 GeV MC
    // m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/22/59/77/18/PDSPProd4a_MC_0.5GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_04.root"); // 0.5 GeV MC
    // m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/54/29/PDSPProd4_data_2GeV_reco2_ntuple_v09_42_03_01.root"); // 2 GeV data
    // m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/52/19/PDSPProd4_data_1GeV_reco2_ntuple_v09_41_00_04.root"); // 1 GeV data
    // m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/AlternateSCE_RITM1506913/00/00/52/35/PDSPProd4_data_1GeV_reco2_ntuple_AltSCEData.root"); // 1 GeV data AltSCE
    // m.AddFile("/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/input/PDSPProd4_data_1GeV_reco2_ntuple_v09_41_00_04.root"); // 1 GeV data local macbook
    // m.AddFile("/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/input/PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"); // 1 GeV MC local macbook
    // m.AddFile("/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/root/PDSPProd4_data_0.5GeV_reco2_ntuple_v09_41_00_04.root"); // 0.5 GeV data local macbook
    // m.AddFile("/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/root/PDSPProd4a_MC_0.5GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_04.root"); // 0.5 GeV MC local macbook
    // m.AddFile("/Users/sungbino/Study/FNAL/ProtoDUNE/samples/PDSPProd4_data_0.5GeV_reco2_ntuple_v09_41_00_04.root"); // 0.5 GeV data local macbook
    // m.AddFile("/Users/sungbino/Study/FNAL/ProtoDUNE/samples/PDSPProd4a_MC_0.5GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_04.root"); // 0.5 GeV MC local macbook
    // m.AddFile("/Users/sungbino/Study/FNAL/ProtoDUNE/samples/pduneana_0.5GeV_030325.root"); // additional 0.5 GeV MC from Jake local macbook
    // m.AddFile("/Users/sungbino/Study/FNAL/ProtoDUNE/samples/pduneana_0.5GeV_031925.root");

    // -- Active 0.5 GeV samples (MF 20260413)
    TString outpath;
    if (sample == "Data") {
        m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/58/25/PDSPProd4_data_0.5GeV_reco2_ntuple_v09_41_00_04.root");
        outpath = "hists_Data_0.5GeV_pionqe0p5.root";
    } else {  // MC
        m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/22/59/77/18/PDSPProd4a_MC_0.5GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_04.root");
        outpath = "hists_MC_0.5GeV_pionqe0p5.root";
        // outpath = "hists_MC_0.5GeV_pionqe0p5_KEff_subt_30MeV.root";
    }

    m.SetOutfilePath(outpath);
    m.Init();
    m.initializeAnalyzer();
    m.initializeAnalyzerTools();
    m.SwitchToTempDir();
    m.Loop();
    m.WriteHist();

    printf("Finished %s: wrote %s\n", sample.Data(), outpath.Data());
}

// ============================================================
// Main entry point. By default runs both MC and Data.
//
// Usage:
//   root -l -b -q run_pionqe.C            // both (default)
//   root -l -b -q 'run_pionqe.C(0)'       // both, explicit
//   root -l -b -q 'run_pionqe.C(1)'       // MC only
//   root -l -b -q 'run_pionqe.C(2)'       // Data only
// ============================================================
void run_pionqe(int which = 0) {
    gSystem->Load("./lib/libDataFormats.so");
    gSystem->Load("./lib/libAnalyzerTools.so");
    gSystem->Load("./lib/libAnalyzers.so");

    bool do_mc = (which == 0 || which == 1);
    bool do_data = (which == 0 || which == 2);
    if (!do_mc && !do_data) {
        printf("ERROR: invalid `which`=%d. Use 0 (both), 1 (MC), or 2 (Data).\n",
               which);
        return;
    }

    if (do_mc) run_one_sample("MC");
    if (do_data) run_one_sample("Data");

    printf("\n========================================================\n");
    printf("  All requested runs complete.\n");
    printf("========================================================\n");
}
