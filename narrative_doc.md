# Narrative Configuration Document

This is going to be a structured, formatted document which describes in narrative detail the elements that go into developing, constructing, and calibrating a model for some application. It is going to be the equivalent of a campaign file, but written out in prose instead of in json code.

The purpose of this document will be twofold:
1. The robot will ingest this document and use it as a guide to building out the actual python code which we need to run the simulation scenarios required for answering a particular scientific question.
2. The document, being a prose record of the elements, features, methodological choices, and data for parameterization and calibration, can be dropped into the methods sections of whatever white paper or academic manuscript which will be the eventual product of the modeling work.


# To Do List

* [ ] Example for training 
  * [ ] Go through EMOD-HIV eswatini configuration file line-by-line, assemble inventory of things that go in
  * [ ] Create narrative document example
* [ ] Read up on plugins and skills, because if the purpose will be to develop AI skill
* [ ] Inventory of skills
  * [ ] Transformation of narrative document into HIVSim
  * [ ] Transformation of EMOD-HIV configuration files into narrative document
    * [ ] Skill which collapses EMOD-HIV configuration files into a python file which generates the json, as compression of thousands of lines of json into something much more compact
  * [ ] Transformation of EMOD-HIV configuration files directly into HIVSim
  * [ ] Need a sanity check: was anything hallucinated along the way?
  * [ ] The robot should be able to generate the requisite HIVSim python code, but may also need to ask questions or make suggestions to the user to make sure the documentation is accurate.


# Narrative document 

## Abstract

The abstract will be structured like a conference abstract, but can be longer and more detailed. It will be the research pitch to the robot. The robot will use this as a guide for building and executing on 

### Background 

This section will contain a prose narrative of the history of the HIV control programs and the healthcare system in the setting that the user is working in. It will include a summary of the scientific question(s) the user wishes to pose.

### Methods/Scenarios

This section will describe the methodology - in particular, it will describe the set of scenarios the user needs to complete the question(s) posed in the background section.

### Results

This section will describe how results are obtained from the simulation - what kinds of postprocessing of the data are required, which metrics should be generated and stored for further analysis. This section can also describe the sorts of figures which the user would like to generate. This section can also describe the sanity checks that the user would conduct to verify that the results are consistent with their intuition.

## Parameterization

This is a list of model parameters which the user would specify. This is also a good place to include scientific papers or other documents for citing these parameters. The user may include the papers or other documents in an accompanying folder  

## Calibration

This is a list of data, preferably cleanly stored in csv tables, which the user would specify as calibration targets. The tables are then stored in a directory which accompany this document

It may be that the user is already using a previously-calibrated model, and they should specify that here.

## Full inventory of documents and directory

* Project Narrative Folder
  * Abstract (markdown document)
  * Directory: Datasets
    * Calibration data
    * Parameterization data
  * Directory: Sources
    * Documents containing supporting references for the project
