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

// Core OpenFOAM includes
#include "fvCFD.H"

// Standard setup includes:
// These define: args, runTime, controlDict, mesh, pimple
#include "setRootCase.H"       // Defines args (if not already), rootCase, caseName
#include "createTime.H"        // Defines runTime, controlDict
#include "createDynamicFvMesh.H" // Defines mesh
#include "createDyMControls.H"   // Defines pimple (via createPimpleControl.H which is included in createControl.H)

// Solver-specific models and utilities
#include "turbulentFluidThermoModel.H" // Should define rhoThermo
#include "radiationModel.H"            // For radiation models if used
#include "pressureControl.H"           // For pressure referencing utilities
#include "fvOptions.H"                 // For fvOptions framework
#include "CorrectPhi.H"                // For phi correction utilities
#include "localEulerDdtScheme.H"       // Specific time scheme
#include "fvcSmooth.H"                 // For smoothing operations

// Custom communication header
#include "comms.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // argList args(argc, argv); // REMOVED: Handled by setRootCase.H
    // if (!args.checkRootCase())
    // {
    //     Foam::FatalError.exit();
    // }

    // Foam::pimpleControl pimple(mesh); // REMOVED: Handled by createDyMControls.H

    // Note for argList is usually added after args is available from setRootCase.H
    // For OpenFOAM.org versions, args might be accessible directly.
    // For OpenFOAM.com versions, it might be through runTime.args().
    // Assuming 'args' is directly available after setRootCase.H for addNote and FatalErrorIn.
    argList::addNote
    (
        "Transient solver for buoyant, turbulent fluid flow"
        " of compressible fluids, including radiation,"
        " with optional mesh motion and mesh topology changes,"
        " and external socket coupling (Forces, Temp, WallHeatFlux)."
    );

    #include "postProcess.H" // For functionObjects, after controlDict and mesh are ready

    // --- Read Solid Properties ---
    Info << "Reading solid properties\n" << endl;
    IOdictionary solidPropsDict
    (
        IOobject
        (
            "solidProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // do not register
        )
    );

    // Hot Block Properties
    const dimensionSet dimSpecificHeat(0, 2, -2, -1, 0, 0, 0); // User's workaround
    const dimensionedScalar hotBlockMass("hotBlockMass", dimMass, solidPropsDict.subDict("hotBlock").lookup("mass"));
    const dimensionedScalar hotBlockCp("hotBlockCp", dimSpecificHeat, solidPropsDict.subDict("hotBlock").lookup("Cp"));
    scalar initialHotBlockTemp = readScalar(solidPropsDict.subDict("hotBlock").lookup("initialTemp"));
    scalar hotBlockTemp = initialHotBlockTemp;

    // Cold Block Properties
    const dimensionedScalar coldBlockMass("coldBlockMass", dimMass, solidPropsDict.subDict("coldBlock").lookup("mass"));
    const dimensionedScalar coldBlockCp("coldBlockCp", dimSpecificHeat, solidPropsDict.subDict("coldBlock").lookup("Cp"));
    scalar initialColdBlockTemp = readScalar(solidPropsDict.subDict("coldBlock").lookup("initialTemp"));
    scalar coldBlockTemp = initialColdBlockTemp;

    // Socket settings
    // controlDict should be available from createTime.H
    word serverIP = controlDict.lookupOrDefault<word>("serverIP", "127.0.0.1");
    int serverPort = controlDict.lookupOrDefault<int>("serverPort", 8080);
    scalar controlInterval = controlDict.lookupOrDefault<scalar>("controlInterval", 1.0);

    wordList wallPatchNames = controlDict.lookupOrDefault<wordList>("wallHeatFluxPatches", wordList());
    if (wallPatchNames.empty()) {
        Info << "No 'wallHeatFluxPatches' specified in controlDict." << endl;
    } else {
        Info << "Will calculate/send wall heat flux for patches: " << wallPatchNames << endl;
    }

    Info << "Attempting to initialize socket connection to "
         << serverIP << ":" << serverPort << endl;
    if (!init_socket(serverIP.c_str(), serverPort))
    {
         Warning << "Socket initialization failed. Continuing without external coupling." << endl;
    }

    const label hotBlockPatchID = mesh.boundary().findPatchID("hotBlockInterface");
    const label coldBlockPatchID = mesh.boundary().findPatchID("coldBlockInterface");

    // Use runTime.args() if 'args' is not directly in scope for FatalErrorIn
    if (hotBlockPatchID < 0) FatalErrorIn(runTime.args().executable()) << "Patch hotBlockInterface not found!" << exit(FatalError);
    if (coldBlockPatchID < 0) FatalErrorIn(runTime.args().executable()) << "Patch coldBlockInterface not found!" << exit(FatalError);

    // Fields are created here. pimple object is now available for pimple.dict() calls.
    // args object is also available for thermo.validate(args.executable(), ...)
    #include "createFields.H" // This file is expected to create T, turbulence, MRF etc.
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    labelList wallPatchIDsToSend;
    if (!wallPatchNames.empty()) {
         wallPatchIDsToSend.setSize(wallPatchNames.size());
         const fvBoundaryMesh& patches = mesh.boundary();
         forAll(wallPatchNames, i) {
             label patchID = patches.findPatchID(wallPatchNames[i]);
             if (patchID < 0) {
                 FatalErrorInFunction // This macro uses the function name
                     << "Cannot find patch named '" << wallPatchNames[i]
                     << "' specified in controlDict entry 'wallHeatFluxPatches'."
                     << exit(FatalError);
             }
             wallPatchIDsToSend[i] = patchID;
         }
    }

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Info<< "\nStarting time loop\n" << endl;
    scalar last_control_time = -controlInterval;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        volScalarField::Boundary& Tbf = T.boundaryFieldRef();

        if (isA<fixedValueFvPatchScalarField>(Tbf[hotBlockPatchID])) {
            fixedValueFvPatchScalarField& hotFixedT =
                dynamic_cast<fixedValueFvPatchScalarField&>(Tbf[hotBlockPatchID]);
            hotFixedT = hotBlockTemp;
        } else {
            FatalErrorInFunction << "hotBlockInterface T BC is not fixedValue!" << exit(FatalError);
        }

        if (isA<fixedValueFvPatchScalarField>(Tbf[coldBlockPatchID])) {
            fixedValueFvPatchScalarField& coldFixedT =
                dynamic_cast<fixedValueFvPatchScalarField&>(Tbf[coldBlockPatchID]);
            coldFixedT = coldBlockTemp;
        } else {
            FatalErrorInFunction << "coldBlockInterface T BC is not fixedValue!" << exit(FatalError);
        }

        autoPtr<volScalarField> divrhoUptr;
        if (correctPhi)
        {
            divrhoUptr.reset
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

        bool performCommunication = false;
        if (runTime.value() > last_control_time + controlInterval - SMALL) {
            if (sock_client && sock_client->isConnected()) {
                 performCommunication = true;
                 Info << "Attempting external coupling communication at t = " << runTime.value() << endl;
            } else {
                 Warning << "Skipping external coupling: Socket not connected." << endl;
                 if (!sock_client || (sock_client && !sock_client->isConnected())) {
                    Info << "Attempting to re-initialize socket connection to "
                         << serverIP << ":" << serverPort << endl;
                    if (!init_socket(serverIP.c_str(), serverPort)) {
                        Warning << "Socket re-initialization failed." << endl;
                    }
                 }
            }
            last_control_time = runTime.value();
        }

        if (performCommunication) {
            auto force_field_data = request_field(runTime.value());
            if (!force_field_data.empty()) {
                if (force_field_data.size() == static_cast<size_t>(F.size())) {
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

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
                if (mesh.changing()) {
                    gh = ((g) & mesh.C()) - ghRef;
                    ghf = ((g) & mesh.Cf()) - ghRef;
                    // Check if MRF object (from createMRF.H in createFields.H) is active
                    if (MRF.active()) MRF.update(); // Changed from MRF.valid()
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

            #include "UEqn.H"
            #include "EEqn.H"

            while (pimple.correct()) {
                #include "pEqn.H"
            }

            if (pimple.turbCorr()) {
                turbulence->correct();
            }
        }

        scalar QnetHot = 0.0;
        scalar QnetCold = 0.0;
        scalar totalPatchAreaHot = 0.0;
        scalar totalPatchAreaCold = 0.0;

        const volScalarField alphaEff(turbulence->alphaEff());
        const volScalarField rho_field(thermo.rho());
        const volScalarField Cp_field(thermo.Cp());

        if (hotBlockPatchID >= 0) {
            const fvPatch& patch = mesh.boundary()[hotBlockPatchID];
            const fvPatchScalarField& alphaEffPatchField = alphaEff.boundaryField()[hotBlockPatchID];
            const fvPatchScalarField& rhoPatchField = rho_field.boundaryField()[hotBlockPatchID];
            const fvPatchScalarField& TPatchField = T.boundaryField()[hotBlockPatchID];
            tmp<scalarField> tSnGradValues = TPatchField.snGrad();
            const scalarField& snGradValues = tSnGradValues();
            forAll(patch.faceCells(), faceI) {
                label faceCell = patch.faceCells()[faceI];
                scalar kEffFace = rhoPatchField[faceI] * Cp_field[faceCell] * alphaEffPatchField[faceI];
                QnetHot += kEffFace * snGradValues[faceI] * patch.magSf()[faceI];
            }
            reduce(QnetHot, sumOp<scalar>());
            totalPatchAreaHot = gSum(patch.magSf());
            if (hotBlockMass.value() * hotBlockCp.value() > VSMALL) {
                hotBlockTemp += (QnetHot * runTime.deltaTValue()) / (hotBlockMass.value() * hotBlockCp.value());
            }
        }

        if (coldBlockPatchID >= 0) {
            const fvPatch& patch = mesh.boundary()[coldBlockPatchID];
            const fvPatchScalarField& alphaEffPatchField = alphaEff.boundaryField()[coldBlockPatchID];
            const fvPatchScalarField& rhoPatchField = rho_field.boundaryField()[coldBlockPatchID];
            const fvPatchScalarField& TPatchField = T.boundaryField()[coldBlockPatchID];
            tmp<scalarField> tSnGradValues = TPatchField.snGrad();
            const scalarField& snGradValues = tSnGradValues();
            forAll(patch.faceCells(), faceI) {
                label faceCell = patch.faceCells()[faceI];
                scalar kEffFace = rhoPatchField[faceI] * Cp_field[faceCell] * alphaEffPatchField[faceI];
                QnetCold += kEffFace * snGradValues[faceI] * patch.magSf()[faceI];
            }
            reduce(QnetCold, sumOp<scalar>());
            totalPatchAreaCold = gSum(patch.magSf());
            if (coldBlockMass.value() * coldBlockCp.value() > VSMALL) {
                coldBlockTemp += (QnetCold * runTime.deltaTValue()) / (coldBlockMass.value() * coldBlockCp.value());
            }
        }

        if (performCommunication) {
            if (!send_temperature(T, runTime.time().value())) {
                Warning << "Failed to send temperature field for t = " << runTime.time().value() << endl;
            }
            std::vector<float> wallFluxDataToSend;
            bool sendHot = false, sendCold = false;
            forAll(wallPatchNames, i) {
                if (wallPatchNames[i] == "hotBlockInterface") sendHot = true;
                if (wallPatchNames[i] == "coldBlockInterface") sendCold = true;
            }
            if(sendHot) {
                scalar avgFluxHot = (totalPatchAreaHot > SMALL) ? QnetHot / totalPatchAreaHot : 0.0;
                wallFluxDataToSend.push_back(static_cast<float>(avgFluxHot));
            }
            if(sendCold) {
                scalar avgFluxCold = (totalPatchAreaCold > SMALL) ? QnetCold / totalPatchAreaCold : 0.0;
                wallFluxDataToSend.push_back(static_cast<float>(avgFluxCold));
            }
            if (!wallFluxDataToSend.empty()) {
                if (!send_wall_heat_flux(wallFluxDataToSend, runTime.time().value())) {
                    Warning << "Failed to send wall heat flux data for t = " << runTime.time().value() << endl;
                }
            }
        }

        rho = thermo.rho();

        if (runTime.write())
        {
            // Info << "Hot Block Temp: " << hotBlockTemp << " K, Cold Block Temp: " << coldBlockTemp << " K" << endl;
            // Info << "QnetHot: " << QnetHot << " W, QnetCold: " << QnetCold << " W" << endl;
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;
    if (sock_client) {
        sock_client->disconnect();
    }
    return 0;
}
// ************************************************************************* //
