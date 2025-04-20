// #############################################################################
// # File: buoyantPimpleFoam.C - Modifications for socket init and temp send #
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
    and mesh topology changes. Includes socket communication for external coupling.

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "comms.H" // Include the updated communications header


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation,"
        " with optional mesh motion and mesh topology changes,"
        " and external socket coupling." // Updated description
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"

    // --- Initialize Socket Connection ---
    // Reads IP and Port from system/controlDict, defaults to 127.0.0.1:8080
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

    word serverIP = controlDict.lookupOrDefault<word>("serverIP", "127.0.0.1");
    int serverPort = controlDict.lookupOrDefault<int>("serverPort", 8080);
    scalar controlInterval = controlDict.lookupOrDefault<scalar>("controlInterval", 5.0); // Use controlDict for interval

    Info << "Attempting to initialize socket connection to "
         << serverIP << ":" << serverPort << endl;
    if (!init_socket(serverIP.c_str(), serverPort))
    {
        // Decide whether to exit or continue without coupling
         Warning << "Socket initialization failed. Continuing without external coupling." << endl;
         // or FatalErrorInFunction << "Socket initialization failed." << exit(FatalError);
    }
    // ----------------------------------

    #include "createFields.H" // F field is created here
    #include "createFieldRefs.H" // psi is referenced here
    #include "initContinuityErrs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    // // Removed patch check - adapt if needed for your specific case
    // if (hotEndPatchID < 0)
    // {
    //     FatalErrorInFunction
    //         << "Cannot find patch named 'hotEnd' in the mesh boundary." << nl
    //         << "Please check boundary file and patch name."
    //         << exit(FatalError);
    // }
    // else
    // {
    //     Info << "Found patch 'hotEnd' with ID: " << hotEndPatchID << endl;
    // }

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // gradT field for writing (if needed)
    // volVectorField gradT = fvc::grad(thermo.T()); // Calculate initial gradient if needed immediately
    // gradT.write(); // Write initial field if needed

    scalar last_control_time = -controlInterval; // Ensure control happens on first step if time > 0

    while (runTime.run())
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh
        // so that it can be mapped and used in correctPhi
        // to ensure the corrected phi has the same divergence
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


        // --- External Coupling Communication ---
        // Check if enough time has passed and if the socket is connected/can connect
        bool sendTempNeeded = false; // Initialize flag for sending temperature
        if (runTime.value() > last_control_time + controlInterval)
        {
             if (sock_client && sock_client->isConnected()) // Check if client exists and is connected
             {
                Info << "Attempting external coupling communication at t = " << runTime.value() << endl;

                // 1. Request Force Field F
                auto force_field_data = request_field(runTime.value());

                if (!force_field_data.empty())
                {
                     if (force_field_data.size() == F.size())
                     {
                        Info << "Received " << force_field_data.size() << " force vectors. Updating F field." << endl;
                        forAll(F, cellI)
                        {
                            // Assuming force_field_data is ordered correctly matching cell IDs
                            F[cellI] = Foam::vector // Use Foam::vector constructor
                            (
                                force_field_data[cellI][0],
                                force_field_data[cellI][1],
                                force_field_data[cellI][2]
                            );
                        }
                        // Optional: Write the updated F field if needed for debugging
                        // F.write();
                     }
                     else
                     {
                         Warning << "Received force field data size (" << force_field_data.size()
                                 << ") does not match internal field size (" << F.size() << "). Skipping update." << endl;
                     }
                }
                else
                {
                    Warning << "Failed to receive valid force field data from external program." << endl;
                    // Decide how to handle this - continue with old F? Zero F? Stop?
                    // For now, it just continues with the existing F field.
                }


                // 2. Send Temperature Field T (after PIMPLE loop converges below)
                // We will send the temperature *after* the PIMPLE loop finishes for this time step.
                // Set a flag or store the time to indicate sending is needed.
                sendTempNeeded = true; // Set flag to true


                last_control_time = runTime.value(); // Update time only if communication was attempted

             } else {
                  Warning << "Skipping external coupling: Socket not connected." << endl;
                  // Attempt to reconnect for the next interval? init_socket() handles this.
             }
        }
        // -------------------------------------


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU.reset(new volVectorField(IOobject::groupName("rhoU", U.group()), rho*U));
                }

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    gh = ((g) & mesh.C()) - ghRef;
                    ghf = ((g) & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        #include "../pimpleFoam/correctPhi.H"

                        // Make the fluxes relative to the mesh-motion
                        fvc::makeRelative(phi, rho, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if (pimple.firstIter() && !pimple.SIMPLErho())
            {
                #include "rhoEqn.H"
            }

            #include "UEqn.H" // F field is used here

            #include "EEqn.H" // Temperature (he) is solved here, thermo.correct() updates T

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H" // Pressure is corrected
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        } // End PIMPLE loop

        // --- Send Temperature Field (if needed after PIMPLE loop) ---
         if (sendTempNeeded && sock_client && sock_client->isConnected())
         {
             Info << "Sending updated temperature field at t = " << runTime.value() << endl;
             if (!send_temperature(thermo.T(), runTime.value())) // Pass the current temperature field T
             {
                 Warning << "Failed to send temperature field for t = " << runTime.value() << endl;
             }
         }
        // ---------------------------------------------------------


        rho = thermo.rho(); // Update rho field based on latest thermo state

        // Optional: Update and write gradT if needed for post-processing
        // gradT = fvc::grad(thermo.T());

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    // Disconnect socket explicitly (though unique_ptr destructor handles it too)
    if (sock_client) {
        sock_client->disconnect();
    }

    return 0;
}


// ************************************************************************* //
