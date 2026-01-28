# Restricted/unrestricted? RI/NORI? RIJCOSX? FROZENCORE? DEFGRIDn(1-3)? Convergence thresholds? DIIS?

TranslationDict = {
    "SP" : "Single Point Energy",
    "OPT" : "Geometry Optimisation",
    "COPT" : "Optimisation in Cartesian Coordinates",
    "ENGRAD" : "Energy and Gradient",
    "NUMGRAD" : "Numerical Gradient",
    "FREQ" : "Frequency",
    "NUMFREQ" : "Numerical Frequency",
    "MD" : "Molecular Dynamics",
    "TD-DFT" : "Time-Dependent Density Functional Theory",
    "OPT+FREQ" : "Geometry Optimisation followed by Frequency",

    "CCSD" : "Coupled-Cluster method with single and double excitations",
    "CCSD(T)" : "Coupled-Cluster method with single, double and perturbative triple excitations",
    "DLPNO-CCSD" : "Domain based Local Pair Natural Orbital Coupled-Cluster method with single and double excitations (closed-shell only)",
    "DLPNO-CCSD(T)" : "Domain based Local Pair Natural Orbital Coupled-Cluster method with single, double and perturbative triple excitations",
    "HF" : "Hartree Fock Theory",
    "DFT" : "Density Functional Theory",
    "MP2" : "Møller–Plesset 2nd order perturbation theory",

    "LDA" : "Local Density Approximation (Not Recommended)",
    
    "BP86" : "Becke '88 exchange and Perdew '86 correlation",

    "B3LYP" : "a GGA hybrid functional (20% HF exchange)",
    "PBE0" : "a GGA hybrid functional",

    "M06L" : "a meta-GGA functional",
    "r2SCAN" : "the regularised and restored SCAN functional (meta-GGA)",
    "B97M-V" : "Head-Gordon's meta-GGA functional with VV10 nonlocal correlation",
    "B97M-D4" : "the modified version of B97M-V (meta-GGA) with DFT-D4 correction",

    "M06" : "a meta-GGA hybrid functional (27% HF exchange)",
    "M062X" : "a meta-GGA hybrid functional (54% HF exchange)",
    "PW6B95" : "a meta-GGA hybrid functional from Truhlar",

    "wB97" : "Head-Gordon’s fully variable range-separated DF 𝜔B97",
    "wB97X" : "Head-Gordon’s range-separated DF 𝜔B97X with minimal Fock exchange",
    "wB97X-D4" : "the modified version of range-separated 𝜔B97X-V with DFT-D4 correction",

    "wB97M-V" : "Head-Gordon’sDF 𝜔B97M-V with VV10 nonlocal correlation",
    "wB97M-D4" : "the modified version of 𝜔B97M-V with DFT-D4 correction",

    "B2PLYP" : "a perturbatively corrected double-hybrid functional - Grimme's mixture of B88, LYP and MP2",
    "wB2PLYP" : "a range-separated double-hybrid functional, with the correlation contributions based on B2PLYP, optimized for excitation energies"
    }

class Orca:
    """
    Data taken from ORCA manual Release 6.0 dated December 13, 2024
    """

    ORCA_RUNTYPES = ["SP", "OPT", "FREQ", "OPT+FREQ", "NUMFREQ", "TD-DFT"]

    ORCA_METHODS = ["DFT", "HF", "MP2", "CCSD", "CCSD(T)", "DLPNO-CCSD","DLPNO-CCSD(T)"]

# Note: Require an "Other" option for the user to input their desired un-listed functional
    ORCA_FUNCTIONALS = [
                    "BP86", "B3LYP", "PBE0", "M06L", "r2SCAN", "B97M-V", "B97M-D4", 
                    "M06", "M062X", "PW6B95", "wB97", "wB97X", "wB97X-D4rev", "wB97X-D4", 
                    "wB97M-V", "wB97M-D4rev", "wB97M-D4", "B2PLYP", "wB2PLYP", "LDA"
                    ]
            
    ORCA_BASIS_SETS = [
                    '6-31G(d)', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ',
                    'aug-cc-pVTZ', 'aug-cc-pVQZ', 'def2-SVP', 'def2-TZVP', 'def2-TZVP(-f)', 
                    'def2-QZVP', 'def2-TVZPP', 'def2-QZVPP', 'def2-TZVPPD', 'def2-QZVPPD', 
                    'ma-def2-SVP', 'ma-def2-TZVP', 'ma-def2-QZVP',
                    ]
    
    ORCA_DISPERSION_CORRECTIONS = ["", "D4rev", "D4", "D3BJ", "D3ZERO"]
    
    ORCA_SOLVENT_MODELS = ["", "CPCM", "SMD", "COSMO-RS", "DRACO"]

    ORCA_SOLVENTS = [
            "", "Water", "Methanol", "Ethanol", "Acetone", "DMSO", 'Toluene', 'Aniline', 
            'Benzene', 'Chloroform', 'Carbon Disulfide', 'DCM', 'diethyl ether', 'DMF', 
            'Ethyl Acetate', 'Nitromethane', 'THF',
        ]

class Gaussian:
    """
    ...
    """

    GAUSSIAN_RUNTYPES = []

    GAUSSIAN_METHODS = ["HF", "DFT", "MP2", "CCSD", "CCSD(T)",]
            
    GAUSSIAN_FUNCTIONALS =[
                        'APFD', 'B3LYP', 'BPV86', 'B3PW91', 'CAM-B3LYP','HCTH', 'HSEH1PBE', 
                        'LSDA', 'MPW1PW91', 'PBEPBE', 'TPSSTPSS', 'WB97XD'
                        ]
            
    GAUSSIAN_BASIS_SETS = [
                        'STO-3G', '3-21G', '6-31G(d)', '6-31G(d,p)', 'LANL2DZ',
                        'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-p5VZ', 'cc-p6VZ',
                        'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-p5VZ', 
                        'aug-cc-p6VZ', 'Def2SV', 'Def2TZV', 'Def2QZV', 'Def2SVP', 
                        'Def2TZVP', 'Def2QZVP', 'Def2SVPP', 'Def2TZVPP', 
                        'Def2QZVPP', 'Def2TZVPD', 'Def2QZVPD',
                        ]
    
    GAUSSIAN_SOLVENT_MODELS = ["None", "PCM", "SMD", "IEFPCM", "CPCM"]

    GAUSSIAN_SOLVENTS = [
            "None", "Water", "Methanol", "Ethanol", "Acetone", "DMSO", "Benzene", "Chloroform", 
            "Toluene", "Acetonitrile", "Dichloromethane", "DiethylEther", "Hexane", "Heptane", 
            "Octanol", "THF", "DMF", "Ethyl Acetate", "Nitromethane"
        ]
    