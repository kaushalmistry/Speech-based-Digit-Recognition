========================================================================
    CONSOLE APPLICATION : 224101031_DigitRecognition Project Overview
========================================================================

AppWizard has created this 224101031_DigitRecognition application for you.

This file contains a summary of what you will find in each of the files that
make up your 224101031_DigitRecognition application.


224101031_DigitRecognition.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard.
    It contains information about the version of Visual C++ that generated the file, and
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

224101031_DigitRecognition.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

224101031_DigitRecognition.cpp
    This is the main application source file.

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named 224101031_DigitRecognition.pch and a precompiled types file named StdAfx.obj.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////


The whole code is the collection of the modular functions
=> There are major 5 functions which are used in the main function and those are listed below:
    -> generateUniverse()
        - The function reads the training files and generates the universe out of it and stores it in the universe.txt file
    
    -> generateCodeBook()
        - Using the LBG algorithm it generates the codebook of size 32 using the universe file generated at above function.

    -> trainHMM()
        - Using the training file (i.e. first 25 utterances of the all digits) it trains the HMM Models and generates the final models for every digits
    
    -> testHMM()
        - It tests all the test files (i.e. 26 to 30 utterances of all the digits) and shows the accuracy along with the winning probability of all the digits
    
    -> liveTesting()
        - This function records the digit live and recognises the spoken digit
    

    -> readCodebook()
        - this reads the final codebook in the codebook array to generate the observation sequences for all the files

=> the file structure is as follows:
    -> "224101031_dataset" folder is having the dataset used for training and testing

    -> "files" folder is having the input to the HMM (initial model lamda) in the "input" directory
        - The "output" directory is having the output of generateUniverse and generateCodeBook functions i.e. universe file and codebook files
        - The "models" directory is having all the intermediate models and the final models for all the digits 