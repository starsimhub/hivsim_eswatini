# Narrative Documentation

This document is a structured, prose description of the elements that go into development, construction, and calibration of a model for a specific application.

# To Do List

* [ ] Reconcile this with STIsim's actual modules: add additional structure to parameterization and calibration sections, so the robot has to do less work to figure out which parameters inform which elements of HIVSim
* [ ] Determine system for versioning - what happens when this document changes and we need to run the model again? Do we create a new branch, or a new commit? Or something else?


# Template

## Abstract

The abstract is structured like a conference abstract, but can be longer and more detailed. It is the research pitch to the robot, which will use it as a guide for building and executing the simulation study.

### Background

A prose narrative of the history of HIV control programs and the healthcare system in the setting the user is working in. This will serve as an overview of assumptions about the ground truth in the model.

### Scientific Question

A summary of the scientific question(s) the user wants to answer.

### Methods/Scenarios

The methodology — in particular, the set of scenarios needed to answer the question(s) posed in the previous section. This will include references to parameters and other assumptions used to conduct this model. 

### Results

How results are obtained from the simulation: what postprocessing of the data is required, which metrics should be generated and stored for further analysis, what figures the user wants produced, and what sanity checks the user would run to confirm the results match their intuition.

## Parameterization

A list of model parameters the user would specify, ideally with citations to the scientific papers or other documents that justify each one. Supporting papers or documents can live in an accompanying `Sources` folder (see the inventory below).

## Calibration

A list of datasets — preferably as clean CSV tables — that the user specifies as calibration targets. If the user is building on a previously-calibrated model rather than calibrating from scratch, that should be specified here too.

## Full inventory of documents and directories

* Project Narrative Folder
  * Abstract (Markdown document)
  * Directory: Datasets
    * Calibration data
    * Parameterization data
  * Directory: Sources
    * Documents containing supporting references for the project
