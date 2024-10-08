/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
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
    chtMultiRegionFoam

Description
    Solver for steady or transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, buoyancy effects, turbulence,
    reactions and radiation modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "pressureReference.H"
#include "hydrostaticInitialisation.H"
#include "Time.H"
#include <iostream>
#include <fstream>
#include <iomanip>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    pimpleMultiRegionControl pimples(fluidRegions,solidRegions);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFluidPressureControls.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    
    //Keff = 1.0;
    std::fstream outfile;
    //std::fstream outfile1;
    std::fstream outfile2;

    outfile .open("keff.dat",std::fstream::out);
    //outfile1.open("TotalPower.dat"   ,std::fstream::out);
    outfile2.open("keffAdjoint.dat",std::fstream::out);

    while (pimples.run(runTime)  ) //
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++; //


        // Optional number of energy correctors
        const int nEcorr = pimples.dict().lookupOrDefault<int>
        (
            "nEcorrectors",
            1
        );

        #include "setMultiRegionDeltaT.H"
        //#include "updateXS.H"

            
        
        //update cross creations
        //#include "updateCrossSections_2D_TWIGL_linear_perturbation.H"

        // --- PIMPLE loop
        while (pimples.loop() )//
        {
            //List<tmp<fvVectorMatrix>> UEqns(fluidRegions.size());

            //for(int Ecorr=0; Ecorr<nEcorr; Ecorr++) //
            TotalFissionPower *= 0;
            forAll(solidRegions, i)
            {
                TotGroupFissionEnergy[i]*=0;
                forAll(phi[i], g)
                {
                    phiErr[i][g]=1.0;
                    TotGroupFissionEnergy[i] += SigmaF[i][g] * phi[i][g]  * Ef[i][g]  ; //
                }
                TotalFissionPower += gSum(fvc::volumeIntegrate(TotGroupFissionEnergy[i])); 
            }

            //Info<<"All Region total volume integrating power: "<<TotalFissionPower<<" W\n"<<endl;
            //if (outfile1.is_open())
            //{
            //    outfile1<<"Time = "<<std::setprecision(10)<<runTime.userTimeName()<<" All Region total volume integrating power: "<<std::setprecision(10)<<TotalFissionPower<<" W"<<std::endl;
            //    outfile1.close();
            //}  
            //else
            //{
            //    outfile1.open("TotalPower.dat",std::fstream::app);
            //    outfile1<<"Time = "<<std::setprecision(10)<<runTime.userTimeName()<<" All Region total volume integrating power: "<<std::setprecision(10)<<TotalFissionPower<<" W"<<std::endl;
            //    outfile1.close();
            //}

            MaxPhiErr = 1.0;//initial err

            //
            while ( MaxPhiErr >= Powertol )
            {            
                
                volFissions *= 0; //
                prevFissions *= 0;
                TotalFissionPower0 *= 0;
                TotalFissionPower1 *= 0;
                
            
                forAll(solidRegions, i)
                {
                    TotGroupFissions0[i]*=0;
                    TotGroupFissionEnergy[i]*=0; 
                    //
                    forAll(phi[i], g) //Q(n)
                        {
                            TotGroupFissions0[i] +=     NuSigmaF[i][g] * phi[i][g] ; //
                            TotGroupFissionEnergy[i] += SigmaF[i][g] * phi[i][g]  * Ef[i][g] ; //
                            phiPre[i][g] = phi[i][g];
            
                            SigmaR[i][g] *=0;
                            SigmaR[i][g] +=SigmaA[i][g];
                            forAll(phi[i], j)
                            {
                                if (g != j)
                                {
                                    SigmaR[i][g] += SigmaS[i][g][j];
                                }
                                else
                                {
                                    SigmaR[i][g] += SigmaS[i][g][j];
                                }  
                            }
                            
                        }
                    prevFissions += gSum(fvc::volumeIntegrate(TotGroupFissions0[i])); 
                    TotalFissionPower0 += gSum(fvc::volumeIntegrate(TotGroupFissionEnergy[i])); 
                      
                }
                
                //scalarIO.writeOpt(TotalFissionPower);
            
                forAll(solidRegions, i) //
                {
                    Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                                      
                    #include "setRegionSolidFields.H"  //for each region the region elds are set using the last timestep
                    #include "solveSolid.H" //Solve the solid region energy equation  
                }
            
            
                //nMesh *= 0; //Get All Mesh Nnmbers
                forAll(solidRegions, i)
                {
                    //nMesh+=phi[i][0].size();
            
                    TotGroupFissions[i]*=0;
                    TotGroupFissionEnergy[i]*=0; 
            
                    forAll(phi[i], g) //Q(n) 
                      {
                
                        TotGroupFissions[i] += NuSigmaF[i][g] * phi[i][g] ;
                        TotGroupFissionEnergy[i] += SigmaF[i][g] * phi[i][g]   * Ef[i][g] ; //
                        
                        MSE[g] += gSum(phiErr[i][g]);  
                      }
            
                    volFissions += gSum(fvc::volumeIntegrate(TotGroupFissions[i])); 
                    TotalFissionPower1 += gSum(fvc::volumeIntegrate(TotGroupFissionEnergy[i])); 
                    
                    Info<<"Region : "<<solidRegions[i].name()<<" Total volume integrating neutron source: "<<volFissions<<" n/(m3*s)"<<endl;  
                                     
                }   
                //
            
                Info<< "Solve the precursor nucleus concenttration equation "<<  endl; 
                forAll(solidRegions, i)
                {
                    forAll(Precursor[i], k) //
                    {
                        delayedNeutronSource[i][k]=TotGroupFissions[i] * precBeta[i][k] / Keff; // Keff
                        fvScalarMatrix precEqn
                        (
                              fvm::ddt(Precursor[i][k])
                            + fvm::Sp(precLambda[i][k],Precursor[i][k])
                            == 
                              delayedNeutronSource[i][k] 
                        );
                        precEqn.relax(); // relax
                        precEqn.solve(); //solve
                    }
                }
            
                Keff0= Keff.value();
                Keff *= mag(volFissions/prevFissions);
            
                //prevFissions=volFissions;//
                
                Info<< "Effective multiplication factor keff = "<< Keff << endl;                
                
            
                //compute phi max error
                //forAll(MSE, g)
                //{
                //    MSE[g] /= nMesh;
                //    RMSE[g] = Foam::pow(MSE[g],0.5);
                //}
            
                MaxPhiErr =  mag(1-TotalFissionPower1/TotalFissionPower0);
                //MaxPhiErr = max(RMSE);
                Info<< "Total Power residual="<< MaxPhiErr << endl;  
                //if (outfile1.is_open())
                //{
                //    outfile1<<"Total Power residual= "<< std::setprecision(10)<< MaxPhiErr << std::endl;  
                //    outfile1.close();
                //}  
                //else
                //{
                //    outfile1.open("PowerResidual.dat",std::fstream::app);
                //    outfile1<<"Total Power residual= "<< std::setprecision(10)<< MaxPhiErr << std::endl;  
                //    outfile1.close();
                //}
            
                Info<< "Effective multiplication factor keff = "<< Keff << endl;
                if (outfile.is_open())
                {
                    outfile<<"Effective multiplication factor keff = "<<std::setprecision(10)<<Keff.value() << std::endl;
                    outfile.close();
                }  
                else
                {
                    outfile.open("keff.dat",std::fstream::app);
                    outfile<<"Effective multiplication factor keff = "<< std::setprecision(10) <<Keff.value() << std::endl;
                    outfile.close();
                }
                
                
            }
            //end forward flux calculation 

            #include "CalculatePinPower.H"

            //start adjoint flux calculation
            while ( MaxAdjointFissionErr >= Powertol )
            {

                volFissionsAdjoint *= 0; //
                prevFissionsAdjoint  *= 0;                

                forAll(solidRegions, i)
                {
                    ChiAdjointFlux[i]*=0;
                    forAll(phiAdjoint[i], j)
                    {
                        ChiAdjointFlux[i] +=  (1-Beta[i]) * ChiPrompt[i][j]  * phiAdjoint[i][j]; //
                        forAll(Precursor[i], k) 
                        {
                            ChiAdjointFlux[i] += precBeta[i][k] * ChiDelayed[i][k][j] * phiAdjoint[i][j];
                        }

                        SigmaR[i][j] *=0;
                        SigmaR[i][j] +=SigmaA[i][j];
                        forAll(phi[i], g)
                        {
                            if (g != j)
                            {
                                SigmaR[i][j] += SigmaS[i][j][g];
                            }
                            else
                            {
                                SigmaR[i][j] += SigmaS[i][j][g];
                            }  
                        }
                    }
                    prevFissionsAdjoint  += gSum(fvc::volumeIntegrate(ChiAdjointFlux[i])); 
                }

                forAll(solidRegions, i) //
                {
                    Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                                      
                    #include "setRegionSolidFields.H"  //for each region the region elds are set using the last timestep
                    #include "solveSolidAdjointFlux.H" //Solve the solid region energy equation  

                    //
                }

                forAll(solidRegions, i)
                {
                    ChiAdjointFlux[i]*=0;
                    forAll(phiAdjoint[i], g)
                    {
                        ChiAdjointFlux[i] +=  (1-Beta[i]) * ChiPrompt[i][g]  * phiAdjoint[i][g]; //Nu1=2.43
                        forAll(Precursor[i], k) 
                        {
                            ChiAdjointFlux[i] += precBeta[i][k] * ChiDelayed[i][k][g] * phiAdjoint[i][g];
                        }
                    }
                    volFissionsAdjoint  += gSum(fvc::volumeIntegrate(ChiAdjointFlux[i])); 
                }

                KeffAdjoint *= mag(volFissionsAdjoint/prevFissionsAdjoint);
                MaxAdjointFissionErr =  mag(1-volFissionsAdjoint/prevFissionsAdjoint);

                Info<< "Maximum Power Error at Time  = "<< MaxAdjointFissionErr << endl; 
                Info<< "Effective multiplication factor keff Adjoint = "<< KeffAdjoint << endl; 
                if (outfile2.is_open())
                {
                    outfile2<<"Effective multiplication factor keff Adjoint = "<<std::setprecision(10)<<KeffAdjoint.value() << std::endl;
                    outfile2.close();
                }  
                else
                {
                    outfile2.open("keffAdjoint.dat",std::fstream::app);
                    outfile2<<"Effective multiplication factor keff Adjoint = "<< std::setprecision(10) <<KeffAdjoint.value() << std::endl;
                    outfile2.close();
                }
                 
            }

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  runTime = " << runTime.userTimeName() << " s"
            << nl << endl;
        
        Info<< "Time = " << runTime.userTimeName() << nl << endl;
    }

    Info<< "End\n" << endl;
    outfile.close();
    //outfile1.close();
    outfile2.close();

    return 0;
}


// ************************************************************************* //
