#!/bin/bash

foamCleanTutorials; touch case.foam; gmshToFoam ./meshing/cyl_fluid.msh; ./buoyantPimpleFoam
