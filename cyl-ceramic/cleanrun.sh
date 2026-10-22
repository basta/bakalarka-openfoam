#!/bin/bash

foamCleanTutorials; touch case.foam; gmshToFoam ./meshing/mesh_new.msh; ./buoyantPimpleFoam
