import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(""),
    InputDirectoryName = cms.string("muonEffs"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("Effciency_ZMuMu_DoubleVoigPlusExpo_RecoIDIso_Uncorrected_2Dbin.root"),

    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),

    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "81", "101", "GeV/c^{2}"),

        #Variables used for passing probe
        pt = cms.vstring("Probe p_{T}", "10", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.4", "2.4", ""),
        dzPV = cms.vstring("dz(PV)", "-99", "99", ""),
        RelPFIso = cms.vstring("PFIso/p_{T}", "0", "99", ""),
        RelPFIso_PFWeighted = cms.vstring("PFIsoWGT/p_{T}", "0", "99", ""),
        RelPFIso_PUPPI = cms.vstring("PFIsoPUPPI/p_{T}", "0", "99", ""),
        RelTrkIso = cms.vstring("TrkIso/p_{T}", "0", "99", "")
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        passingGlb = cms.vstring("isGLBmuon", "dummy[true=1,false=0]"),
        passingtrkiso = cms.vstring("passingtrkiso", "dummy[true=1, false=0]"),
        passingTight = cms.vstring("passingTight", "dummy[true=1, false=0]"),
        passingHLT_Mu17Mu8dzIsoVVL = cms.vstring("passingHLT_Mu17Mu8dzIsoVVL", "dummy[true=1, false=0]")
    ),

    Cuts = cms.PSet(
        pt20 = cms.vstring("P_{T} > 20GeV", "pt", "20"),
        eta2p4_upper = cms.vstring("#eta < 2.4", "eta", "2.4"),
        eta2p4_lower = cms.vstring("#eta > -2.4", "eta", "-2.4"),
        dzPV_upper = cms.vstring("dz(PV) < 0.5", "dzPV", "0.5"),
        dzPV_lower = cms.vstring("dz(PV) > -0.5", "dzPV", "-0.5"),
        RelPFIso12 = cms.vstring("RelPFIso_Tight", "RelPFIso", "0.12"),
        RelPFIso_WGT12 = cms.vstring("RelPFIso_PFWeighted_Tight", "RelPFIso_PFWeighted", "0.12"),
        RelPFIso_PUPPI12 = cms.vstring("RelPFIso_PUPPI_Tight", "RelPFIso_PUPPI", "0.12"),
        RelTrkIso10 = cms.vstring("RelTrkIso_Tight", "RelTrkIso", "0.10")
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        DoubleVoigPlusExpo = cms.vstring(
            #signals
            "Voigtian::signal1(mass, mean1[9.10297e+01, 86, 96], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[9.10297e+01, 86, 96], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            
            #backgrounds
            "Exponential::backgroundPass(mass, lp[-0.25,-2,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.25,-2,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        IDTight_PFIso = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("passingHLT_Mu17Mu8dzIsoVVL", "true", "dzPV_lower", "above", "dzPV_upper", "below", "passingTight", "true", "RelPFIso12", "below"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(10, 17, 30, 40, 50, 70, 250, 1000),
                eta = cms.vdouble(-2.4, -2.1, -1.9, -1.5, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 1.9, 2.1, 2.4),
            ),
            BinToPDFmap = cms.vstring("DoubleVoigPlusExpo")
		),
    )
)

process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)

