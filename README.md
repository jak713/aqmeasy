## AQMEasy is a Graphical User Interface built for the [AQME (Automated Quantum Mechanical Environments)](https://github.com/jvalegre/aqme) software.

## Installation
1. Clone this repository: `git clone https://github.com/jak713/aqmeasy.git`
2. Create new conda environment: `conda create -n aqmeasy python=3.12`
    Note, if coming from AQME, you can use the same environment already used for AQME with additional dependencies found in `requirements.txt`
3. `conda activate aqmeasy`
4. Install pip dependencies: `pip install -r requirements.txt`
5. Run the program with `python3 main.py`

### TODO:
## CSEARCH
##### - Add the complex_type option when TM is found !
##### - Fix run aqme command  !!!!
##### - expand the pubchem search  ... !
##### ---- fix the warning where the smiles entered are not valid (rn nothing really happens but the change is only regirestered when the smiles are valid (or thats how it should be)), similarly changing smiles removes the view of the constraints but doesnt remove them completely at the index level, which I think should be done?
## QPREP
#### - charge/mult needs to be incorporated into model and passed to aqme properly
#### - look for bugs!
## QCORR
#### - robust file list, with drag n drop
#### - connection with qprep param qwidget for qm_input
#### - robust models design... need to think about this
#### - json parsing from cclib and everything to do with that
## QDESP
#### --- yet to begin