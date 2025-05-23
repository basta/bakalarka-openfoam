// #############################################################################
// # File: buoyantPimpleFoam.C - Added wall heat flux calculation & sending  #
// #############################################################################
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of compressible fluids
    for ventilation and heat-transfer, with optional mesh motion
    and mesh topology changes. Includes socket communication for external coupling,
    sending temperature field and calculated wall heat flux.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
// Removed incorrect include: #include "basicThermophysicalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "comms.H" // Include the updated communications header


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation,"
        " with optional mesh motion and mesh topology changes,"
        " and external socket coupling (Forces, Temp, WallHeatFlux)." // Updated description
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"

    // --- Read Solver Configuration ---
    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Socket settings
    word serverIP = controlDict.lookupOrDefault<word>("serverIP", "127.0.0.1");
    int serverPort = controlDict.lookupOrDefault<int>("serverPort", 8080);
    scalar controlInterval = controlDict.lookupOrDefault<scalar>("controlInterval", 1.0); // Default to 1s

    // Wall patches for heat flux calculation
    wordList wallPatchNames = controlDict.lookupOrDefault<wordList>("wallHeatFluxPatches", wordList());
    if (wallPatchNames.empty()) {
        Info << "No 'wallHeatFluxPatches' specified in controlDict. Wall heat flux will not be calculated or sent." << endl;
    } else {
        Info << "Will calculate and send average wall heat flux for patches: " << wallPatchNames << endl;
    }
    // ----------------------------------


    // --- Initialize Socket Connection ---
    Info << "Attempting to initialize socket connection to "
         << serverIP << ":" << serverPort << endl;
    if (!init_socket(serverIP.c_str(), serverPort))
    {
         Warning << "Socket initialization failed. Continuing without external coupling." << endl;
    }
    // ----------------------------------

    #include "createFields.H" // F field is created here
    #include "createFieldRefs.H" // psi is referenced here
    #include "initContinuityErrs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate(); // Ensure turbulence model is set up

    // --- Get Wall Patch IDs ---
    labelList wallPatchIDs;
    if (!wallPatchNames.empty()) {
         wallPatchIDs.setSize(wallPatchNames.size());
         // FIX: Use fvBoundaryMesh type
         const fvBoundaryMesh& patches = mesh.boundary();
         forAll(wallPatchNames, i) {
             label patchID = patches.findPatchID(wallPatchNames[i]);
             if (patchID < 0) {
                 FatalErrorInFunction
                     << "Cannot find patch named '" << wallPatchNames[i]
                     << "' specified in controlDict entry 'wallHeatFluxPatches'."
                     << exit(FatalError);
             }
             wallPatchIDs[i] = patchID;
         }
    }
    // --------------------------


    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    scalar last_control_time = -controlInterval; // Ensure control happens on first valid step

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU.reset
            (
                new volScalarField
                (
                    IOobject::groupName("divrhoU", runTime.timeName()),
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }

        if (LTS)
        {
            #include "../pimpleFoam/setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        // --- External Coupling Communication Check ---
        bool performCommunication = false;
        if (runTime.value() > last_control_time + controlInterval) {
            if (sock_client && sock_client->isConnected()) {
                 performCommunication = true;
                 Info << "Attempting external coupling communication at t = " << runTime.value() << endl;
            } else {
                 Warning << "Skipping external coupling: Socket not connected." << endl;
                 // Attempt to reconnect for the next interval? init_socket() handles this.
            }
             // Update time only if communication was attempted (successful or not)
             last_control_time = runTime.value();
        }
        // ------------------------------------------


        // --- Request Force Field F (if communicating) ---
        if (performCommunication) {
            auto force_field_data = request_field(runTime.value());
            if (!force_field_data.empty()) {
                 if (force_field_data.size() == static_cast<size_t>(F.size())) {
                    // Info << "Received " << force_field_data.size() << " force vectors. Updating F field." << endl; // Less verbose
                    forAll(F, cellI) {
                        F[cellI] = Foam::vector(
                            force_field_data[cellI][0],
                            force_field_data[cellI][1],
                            force_field_data[cellI][2]
                        );
                    }
                 } else {
                     Warning << "Received force field data size (" << force_field_data.size()
                             << ") does not match internal field size (" << F.size() << "). Skipping update." << endl;
                 }
            } else {
                Warning << "Failed to receive valid force field data from external program." << endl;
            }
        }
        // ---------------------------------------------


        // --- Pressure-velocity PIMPLE corrector loop ---
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid()) {
                    rhoU.reset(new volVectorField(IOobject::groupName("rhoU", U.group()), rho*U));
                }
                mesh.update();
                if (mesh.changing()) {
                    gh = ((g) & mesh.C()) - ghRef;
                    ghf = ((g) & mesh.Cf()) - ghRef;
                    MRF.update();
                    if (correctPhi) {
                        phi = mesh.Sf() & rhoUf();
                        #include "../pimpleFoam/correctPhi.H"
                        fvc::makeRelative(phi, rho, U);
                    }
                    if (checkMeshCourantNo) {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if (pimple.firstIter() && !pimple.SIMPLErho()) {
                #include "rhoEqn.H"
            }

            #include "UEqn.H" // F field is used here
            #include "EEqn.H" // Temperature (he) is solved here, thermo.correct() updates T

            while (pimple.correct()) { // Pressure corrector loop
                #include "pEqn.H"
            }

            if (pimple.turbCorr()) {
                turbulence->correct(); // Updates turbulence fields like alphat, nut
            }
        } // End PIMPLE loop


        // --- Send Data (if communicating) ---
        if (performCommunication) {
             // 1. Send Temperature Field T
             Info << "Sending updated temperature field at t = " << runTime.value() << endl;
             if (!send_temperature(thermo.T(), runTime.value())) {
                 Warning << "Failed to send temperature field for t = " << runTime.value() << endl;
             }

             // 2. Calculate and Send Wall Heat Flux
             if (!wallPatchIDs.empty()) {
                 Info << "Calculating and sending wall heat flux for specified patches..." << endl;
                 std::vector<float> wallFluxData;
                 wallFluxData.reserve(wallPatchIDs.size());

                 // Get necessary fields
                 const volScalarField& T = thermo.T();
                 // FIX: Avoid dangling reference by creating a copy
                 const volScalarField alphaEff = turbulence->alphaEff();
                 const fvBoundaryMesh& boundary = mesh.boundary();
                 // Get reference to the thermo object to access Cp
                 const basicThermo& basicThermo = thermo;
                 // FIX: Get Cp field once
                 const volScalarField Cp = thermo.Cp();


                 scalar totalAreaSum = 0;    // For sanity check
                 scalar totalFluxSum = 0;    // For sanity check

                 forAll(wallPatchIDs, i) {
                     label patchID = wallPatchIDs[i];
                     const fvPatch& patch = boundary[patchID];
                     const scalarField& patchAreas = patch.magSf(); // Face areas
                     // Access boundary fields correctly
                     const fvPatchScalarField& alphaEffPatchField = alphaEff.boundaryField()[patchID];
                     const fvPatchScalarField& TPatchField = T.boundaryField()[patchID];
                     const fvPatchScalarField& rhoPatchField = rho.boundaryField()[patchID];


                     // Calculate surface normal gradient for temperature on the patch
                     // FIX: Use the temporary result directly in the loop below
                     // tmp<surfaceScalarField> tsnGradT = TPatchField.snGrad(); // This returns tmp<scalarField>
                     // const surfaceScalarField& snGradT = tsnGradT();

                     scalar totalPatchFlux = 0.0;
                     scalar totalPatchArea = 0.0;

                     // Heat flux q = -k_eff * grad(T)_n = - (rho*Cp*alphaEff) * snGrad(T)
                     // Calculate k_eff = rho * Cp * alphaEff on the patch faces

                     // Need Cp on patch faces. (Using adjacent cell value from Cp field)
                     // Need p on patch faces (or adjacent cells) for Cp calculation - Not needed if using Cp field directly
                     // FIX: Remove unused pPatchField declaration
                     // const fvPatchScalarField& pPatchField = p.boundaryField()[patchID];


                     // Get the raw snGrad values for the patch faces
                     tmp<scalarField> tSnGradValues = TPatchField.snGrad();
                     const scalarField& snGradValues = tSnGradValues();


                     forAll(patch.faceCells(), faceI) {
                         // Ensure direct access to boundary field values using []
                         scalar rhoFace = rhoPatchField[faceI];
                         scalar alphaEffFace = alphaEffPatchField[faceI];

                         label faceCell = patch.faceCells()[faceI];
                         scalar CpFace = Cp[faceCell]; // Access Cp field value for the cell

                         scalar kEffFace = rhoFace * CpFace * alphaEffFace;
                         scalar qFace = -kEffFace * snGradValues[faceI]; // Flux density (W/m^2)
                         totalPatchFlux += qFace * patchAreas[faceI]; // Integrate flux (W)
                         totalPatchArea += patchAreas[faceI];         // Integrate area (m^2)
                     }

                     scalar avgFluxDensity = (totalPatchArea > SMALL) ? (totalPatchFlux / totalPatchArea) : 0.0;
                     wallFluxData.push_back(static_cast<float>(avgFluxDensity));

                     totalAreaSum += totalPatchArea;
                     totalFluxSum += totalPatchFlux;

                     // Info << "  Patch '" << patch.name() << "' (ID " << patchID << "): Avg Flux = " << avgFluxDensity << " W/m^2" << endl;
                 }
                  Info << "  Calculated Avg Fluxes (W/m^2): " << wallFluxData << endl;
                  // Info << "  Total Flux Sum (W): " << totalFluxSum << ", Total Area Sum (m^2): " << totalAreaSum << endl;


                 // Send the calculated average flux data
                 if (!send_wall_heat_flux(wallFluxData, runTime.value())) {
                     Warning << "Failed to send wall heat flux data for t = " << runTime.value() << endl;
                 }
             }
        }
        // ------------------------------------


        rho = thermo.rho(); // Update rho field based on latest thermo state

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    // Disconnect socket explicitly
    if (sock_client) {
        sock_client->disconnect();
    }

    return 0;
}


// ************************************************************************* //
