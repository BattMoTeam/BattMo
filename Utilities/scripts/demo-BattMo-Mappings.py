# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 22:47:22 2022

@author: simonc
"""

from pathlib import Path
import os
import dlite
from dlite.triplestore import (
    en, Literal, Triplestore,
    EMMO, OWL, RDF, RDFS, SKOS, XSD,
)
from dlite.mappings import instantiate
from dlite.mappings import mapping_route

import pint
from scipy import integrate

import matplotlib.pyplot as plt


# set directory information
thisdir = Path(__file__).resolve().parent
entitydir = thisdir.parent / 'entities'
outputdir = thisdir / 'output'
ontodir = thisdir.parent / 'ontologies'
datadir = thisdir.parent / 'data'

# define DLite data models

datamodel = "BattMoLiIonBatteryDataModel.json"
datamodel_path = os.path.join(entitydir, datamodel)
BattMoLiIonBatteryDataModel = dlite.Instance.from_url(f'json://{datamodel_path}')

datamodel = "BattMoLiIonElectrodeDataModel.json"
datamodel_path = os.path.join(entitydir, datamodel)
BattMoLiIonElectrodeDataModel = dlite.Instance.from_url(f'json://{datamodel_path}')

datamodel = "BattMoLiIonActiveMaterialDataModel.json"
datamodel_path = os.path.join(entitydir, datamodel)
BattMoLiIonActiveMaterialDataModel = dlite.Instance.from_url(f'json://{datamodel_path}')

collection = dlite.Collection()

# create instances
cell    = BattMoLiIonBatteryDataModel(dims=[1])
pe      = BattMoLiIonElectrodeDataModel(dims=[1])
ne      = BattMoLiIonElectrodeDataModel(dims=[1])
nmc     = BattMoLiIonActiveMaterialDataModel(dims=[1])
gr      = BattMoLiIonActiveMaterialDataModel(dims=[1])

# create dlite collection to group instances
collection.add(label='cell', inst=cell)
collection.add(label='positive_electrode', inst=pe)
collection.add(label='positive_electrode_active_material', inst=nmc)
collection.add(label='negative_electrode', inst=ne)
collection.add(label='negative_electrode_active_material', inst=gr)

# define instance hierarchy
cell.positive_electrode = [pe]
cell.negative_electrode = [pe]

cell.positive_electrode[0].active_material = [nmc]
cell.negative_electrode[0].active_material = [gr]

print(cell)

# create the triplestore
ts = Triplestore("rdflib")
ts.parse(f"{ontodir}/battinfo-merged.ttl")

# BattINFO namespace
BATTINFO = ts.bind(
    'battinfo', 'https://big-map.github.io/BattINFO/ontology/BattINFO#')

# Dict mapping prefLabel to IRI
d = {o.value: s for s, o in ts.subject_objects(SKOS.prefLabel)}

# map model quantities
ts.add_mapsTo(d['ElectrodeContinuumModel'], cell, 'positive_electrode')
ts.add_mapsTo(d['ElectrodeContinuumModel'], cell, 'negative_electrode')
ts.add_mapsTo(d['EndOfDischargeVoltage'], cell, 'Ucut')
ts.add_mapsTo(d['ThermodynamicTemperature'], cell, 'initT')
#ts.add_mapsTo(d['StateOfCharge'], cell, 'SOC')

ts.add_mapsTo(d['ActiveElectrochemicalMaterialContinuumModel'], pe, 'active_material')
ts.add_mapsTo(d['CurrentCollectorContinuumModel'], pe, 'current_collector')

ts.add_mapsTo(d['ActiveElectrochemicalMaterialContinuumModel'], ne, 'active_material')
ts.add_mapsTo(d['CurrentCollectorContinuumModel'], ne, 'current_collector')

#ts.add_mapsTo(d['BruggemanCoefficient'], nmc, 'bruggeman_coefficient')
ts.add_mapsTo(d['ThermalConductivity'], nmc, 'thermal_conductivity')
ts.add_mapsTo(d['ElectronicConductivity'], nmc, 'electrical_conductivity')
ts.add_mapsTo(d['SpecificHeatCapacity'], nmc, 'heat_capacity')
ts.add_mapsTo(d['DiffusionCoefficient'], nmc, 'inter_diffusion_coefficient')

#ts.add_mapsTo(d['BruggemanCoefficient'], gr, 'bruggeman_coefficient')
ts.add_mapsTo(d['ThermalConductivity'], gr, 'thermal_conductivity')
ts.add_mapsTo(d['ElectronicConductivity'], gr, 'electrical_conductivity')
ts.add_mapsTo(d['SpecificHeatCapacity'], gr, 'heat_capacity')
ts.add_mapsTo(d['DiffusionCoefficient'], gr, 'inter_diffusion_coefficient')

# query the triplestore
query_text = """
PREFIX map: <http://emmo.info/domain-mappings#>
PREFIX emmo: <http://emmo.info/emmo#>
PREFIX battinfo: <https://big-map.github.io/BattINFO/ontology/BattINFO#>

SELECT *
WHERE {
   ?subject map:mapsTo battinfo:EMMO_b72eb3ad_8935_4420_a64e_6218de31c0d2 .
}
"""

query_result = ts.query(query_text)

for row in query_result:
    print(f"{row.subject}")

# save the collection
collection.save('json', f'{thisdir}/output/battmo_collection.json', 'mode=w')

# serialize the triplestore
ts.serialize(f'{thisdir}/output/mapping_triples.ttl')