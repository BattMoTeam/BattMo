=============
Features
=============

This section explains some convenient features that the application includes.


Use default materials
=====================
The application has some default materials that you can use in your input parameter setup.
These materials are datasets formed from literature. An overview of which materials are available
can be found on the page `Materials and models <http://app.batterymodel.com/Materials_and_models>`_ together with the reference and dataset details.

Define your own materials
=========================
Did you characterize your own material in the lab, found one in literature, or just simply would like to see how changing material characteristics change the
simulation results? Then you can do this by defining your own material. You can do this on the `Simulation <http://app.batterymodel.com/Simulation>`_ page by selecting 'User defined' in the 
a material selectbox. When 'User defined' is selected, an expander will appear where you can fill in your own parameter values. You can still use parameter values from the default materials by copying
the values from the `Materials and models <http://app.batterymodel.com/Materials_and_models>`_ page.

Visualize your geometry
=======================
We've included a feature on the `Simulation <http://app.batterymodel.com/Simulation>`_ page that visualizes the battery cell geometry depending on the component 
thicknesses and porosities, the length, and the width that are defined in the parameter inputs. 

Download your input data
========================
After defining your parameters on the `Simulation <http://app.batterymodel.com/Simulation>`_ page, you can download them in two different formats.
The first format is the JSON Linked Data format in which the data is structured, together with its metadata, according to the `FAIR data principles <https://www.go-fair.org/fair-principles/>`_ and 
the `5 star open data guidelines <https://5stardata.info/en/>`_ in order to improve interoperability.
In the second format, the input data is structured according to the |battmo| format.

Visualize and download your results
===================================
After a simulation has finished, the results can immediatly be Visualized on the `Results <http://app.batterymodel.com/Results>`_ page.
You can choose from several results parameters and visualize them in line or color plots. You can also download you results as an HDF5 file.