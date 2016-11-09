# VIRUS

A variety of scripts that are useful for VIRUS IFU spectroscopic reductions and inspections.

### Authors

* Greg Zeimann, UT Austin
* Karl Gebhardt, UT Austin
* Daniel Farrow, MPE Garching

### Description

#### bspline.py

* johntfoster's python implementation of Bspline basis function via Cox - de Boor algorithm.  Includes memoization for rapid recalculations.

#### calibration_script.py

* reduces calibration frames to produce masterbias, masterdark, mastertrace_twi, masterarc frames, as well as distortion and fiber model files.

#### call_cure.py

* python interface with many CURE functions

#### examine_darks.py

* In progress script for looking at dark frames

#### gain_readnoise.py

* Measure gain and readnoise from bias/flat frames.  Can be fiber flats or full flats.

#### ifu_normalization

* Script to normalizes IFU to IFU through fiberextract and averaging the twighlight frames as a function of wavelength.  Can apply normalization to FE frames.

#### science_scripty.py

* Science reducitons to be run after calibraiton_script.py

#### utils.py

* Useful mathematical functions including biweight_location (with axis argument) 

