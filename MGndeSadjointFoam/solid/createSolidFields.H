// Initialise solid field pointer lists 
#include "nuclearProperties.H" // read data
//PtrList<solidThermo> thermoSolid(solidRegions.size()); //
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size()); //
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size()); //

PtrList<PtrList<volScalarField>> phi(solidRegions.size()); //
PtrList<PtrList<volScalarField>> phiPre(solidRegions.size()); //prestore
PtrList<PtrList<volScalarField>> phiErr(solidRegions.size()); //error

PtrList<PtrList<volScalarField>> Dg(solidRegions.size()); //
PtrList<PtrList<volScalarField>> InverseVel(solidRegions.size()); //inverse of velocity
PtrList<PtrList<volScalarField>> SigmaA(solidRegions.size()); //absorption
PtrList<PtrList<volScalarField>> SigmaR(solidRegions.size()); //remove
PtrList<PtrList<volScalarField>> SigmaF(solidRegions.size()); //fisssion
PtrList<PtrList<PtrList<volScalarField>>> SigmaS(solidRegions.size());//scatter matrix
PtrList<PtrList<volScalarField>> Qs(solidRegions.size()); //scatter term
PtrList<PtrList<volScalarField>> Qf(solidRegions.size()); //fission term

PtrList<PtrList<volScalarField>> ChiPrompt(solidRegions.size()); // chi prompt
PtrList<PtrList<PtrList<volScalarField>>> ChiDelayed(solidRegions.size()); // chi delayed
PtrList<PtrList<volScalarField>> NuSigmaF(solidRegions.size());  // Nu sigmaF
PtrList<PtrList<volScalarField>> NeutronSource(solidRegions.size()); //neutron source
PtrList<PtrList<volScalarField>> ScatteringTerm(solidRegions.size()); //

PtrList<PtrList<volScalarField>> delayedNeutronSource(solidRegions.size()); //
PtrList<PtrList<volScalarField>> GdelayedNeutronSource(solidRegions.size()); //
PtrList<PtrList<volScalarField>> Precursor(solidRegions.size()); //precursor
PtrList<PtrList<volScalarField>> precLambda(solidRegions.size());//lambda
PtrList<PtrList<volScalarField>> precBeta(solidRegions.size());//beta

PtrList<volScalarField> zeroField(solidRegions.size()); //

PtrList<volScalarField> Beta(solidRegions.size()); //total beta
PtrList<PtrList<volScalarField>> Nu(solidRegions.size()); //nu

PtrList<volScalarField> NuRfTot(solidRegions.size()); //
PtrList<PtrList<volScalarField>> Ef(solidRegions.size()); //fission energy

PtrList<volScalarField> TotGroupFissions(solidRegions.size()); //
PtrList<volScalarField> TotGroupFissions0(solidRegions.size()); //

PtrList<volScalarField> TotGroupFissionEnergy(solidRegions.size()); //

dimensionedScalar Keff("Keff", dimless, scalar(KeffIni));
dimensionedScalar MiNiNum("", dimless/dimArea/dimTime, 1e-16);

dimensionedScalar KeffAdjoint("KeffAdjoint", dimless, scalar(KeffIni));
PtrList<volScalarField> ChiAdjointFlux(solidRegions.size()); 
PtrList<PtrList<volScalarField>> phiAdjoint(solidRegions.size()); 
PtrList<PtrList<volScalarField>> ScatteringTermAdjoint(solidRegions.size()); 
PtrList<PtrList<volScalarField>> NeutronSourceAdjoint(solidRegions.size()); 
//const IOobject scalarIO
//(
//    "TotalFissionEnergy",
//    runTime.timeName(),
//    runTime,
//    IOobject::NO_READ,
//    IOobject::AUTO_WRITE    
//
//);

PtrList<dimensionedScalar> MaxPhi(energyGroups);
PtrList<dimensionedScalar> MaxPhi0(energyGroups);
forAll(MaxPhi,i)
{
    word Iname=name(i+1);
    MaxPhi.set(i, new dimensionedScalar("MaxPhi" + Iname , dimless/dimArea/dimTime, 0));
    MaxPhi0.set(i, new dimensionedScalar("MaxPhi" + Iname , dimless/dimArea/dimTime, 0));
}

scalar Keff0;
List<scalar> MSE(energyGroups);
List<scalar> RMSE(energyGroups);
forAll(RMSE,i)
{
    RMSE[i] = 0;
    MSE[i] = 0;
}
label nMesh=0;

scalar volFissions=1.0;
scalar prevFissions=1.0;

scalar volFissionsAdjoint=1.0;
scalar prevFissionsAdjoint=1.0;

scalar TotalFissionPower=1.0;
scalar TotalFissionPower0=1.0;
scalar TotalFissionPower1=1.0;


scalar Kerr=KeffIni;
scalar MaxPhiErr=1.0;
scalar MaxAdjointPhiErr=1.0;
scalar MaxAdjointFissionErr=1.0;

//List for storage Max phi error
//PtrList<PtrList<scalar>> AllPhiMaxErr(solidRegions.size());
//forAll(solidRegions, i)
//{
//    AllPhiMaxErr.set(i,new PtrList<scalar>(energyGroups));
//    forAll(AllPhiMaxErr[i],j)
//    {
//        AllPhiMaxErr[i][j]=1.0;
//    }
//}

// Populate solid field pointer lists
forAll(solidRegions, i) // the i-th region
{

    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    //Info<< "    Adding to thermoSolid\n" << endl;
    //thermoSolid.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );


    Info<< "    Adding zeroField for zero in source item\n" << endl;
    zeroField.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"zeroField",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::NO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("", dimless, 0)
    	)
    );

    Info<< "Region "<< i<< " Adding total group fission neutron\n" << endl;
    TotGroupFissions.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"TotGroupFissions",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::AUTO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("",dimless/dimVolume/dimTime, 0.0)
    	)
    );

    Info<< "Region "<< i<< " Adding total group fission neutron 0\n" << endl;
    TotGroupFissions0.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"TotGroupFissions0",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::AUTO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("",dimless/dimVolume/dimTime, 0.0)
    	)
    );

    Info<< "Region "<< i<< " Adding total group fission energy \n" << endl;
    TotGroupFissionEnergy.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"TotGroupFissionEnergy",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::AUTO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("",dimless/dimVolume/dimTime, 0.0)
    	)
    );

    Info<< "Region "<< i<< " Adding total group fission energy \n" << endl;
    ChiAdjointFlux.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"ChiAdjointFlux",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::AUTO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("",dimless/dimArea/dimTime, 0.0)
    	)
    );




    Info<< "Region "<< i<< " Adding total precursor fraction n" << endl;
    Beta.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"Beta" ,
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::NO_WRITE
    		),
    	solidRegions[i],
        Beta_ref[i]
    	)
    );

    Info<< "    Adding total reactron rate\n" << endl;
    NuRfTot.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"NuRfTot",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::AUTO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("", dimless, 0)
    	)
    );

    Info<< "Region "<< i<< " Adding Nu\n" << endl;
    Nu.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(Nu[i],j)
    {
        Nu[i].set
        (
        	j,
        	new volScalarField
        	(
        		IOobject
        		(
        			"Nu"+Foam::name(j+1),
        			runTime.timeName(),
        			solidRegions[i],
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
        		),
        	solidRegions[i],
            Nu_ref[i][j]
        	//dimensionedScalar("", dimless, 0)
        	)
        );
    }    

    

    Info<< "Region "<< i<< " Adding flux \n" << endl;
    phi.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(phi[i], j) //
        {
            phi[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"phi"+Foam::name(j+1),
                        	runTime.timeName(),
                        	solidRegions[i],
                        	IOobject::MUST_READ, //
                        	IOobject::AUTO_WRITE 
                        ),
                        solidRegions[i]
                        //dimensionedScalar( dimless/dimArea/dimTime)
                )
            );
        }

    Info<< "Region "<< i<< " Adding flux \n" << endl;
    phiAdjoint.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(phiAdjoint[i], j) //
        {
            phiAdjoint[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"phiAdjoint"+Foam::name(j+1),
                        	runTime.timeName(),
                        	solidRegions[i],
                        	IOobject::MUST_READ, //
                        	IOobject::AUTO_WRITE 
                        ),
                        solidRegions[i]
                        //dimensionedScalar( dimless/dimArea/dimTime)
                )
            );
        }


    Info<< "Region "<< i<< " Adding flux Previous \n" << endl;
    phiPre.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(phiPre[i], j) //
        {
            phiPre[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"phiPre"+Foam::name(j+1),
                        	runTime.timeName(),
                        	solidRegions[i],
                        	IOobject::MUST_READ, //
                        	IOobject::AUTO_WRITE 
                        ),
                        solidRegions[i],
                        dimensionedScalar("",dimless/dimArea/dimTime,0)
                )
            );
        }

    Info<< "Region "<< i<< " Adding flux Error \n" << endl;
    phiErr.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(phiErr[i], j) //
        {
            phiErr[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"phiErr"+Foam::name(j+1),
                        	runTime.timeName(),
                        	solidRegions[i],
                        	IOobject::MUST_READ, //
                        	IOobject::AUTO_WRITE 
                        ),
                        solidRegions[i],
                        dimensionedScalar("",dimless,1.0)
                )
            );
        }

    Info<< " Adding flux diffusion coefficient\n" << endl;
    Dg.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(Dg[i], j)
        {
            Dg[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"Dg"+Foam::name(j+1),
                        	runTime.timeName(),
                        	solidRegions[i],
                        	IOobject::NO_READ,
                        	IOobject::NO_WRITE
                        ),
                        solidRegions[i],
                        D_ref[i][j]
                )
            );
        }


    Info<< "Region "<< i<< " Adding velocity\n" << endl;
    InverseVel.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(InverseVel[i], j)
    {
        InverseVel[i].set(
            j,
            new volScalarField(
                IOobject
                    (
                    	"InverseVel"+Foam::name(j+1),
                    	runTime.timeName(),
                    	solidRegions[i],
                    	IOobject::NO_READ,
                    	IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    InverseVelocity_ref[i][j]
            )
        );
    }

    Info<< "Region "<< i<< " Adding absorption cross section \n" << endl;
    SigmaA.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(SigmaA[i], j)
        {
         SigmaA[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"SigmaA"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                     solidRegions[i],
                     SigmaA_ref[i][j]
             )
         );
        }

    Info<< "Region "<< i<< " Adding total loss cross section \n" << endl;
    SigmaR.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(SigmaR[i], j)
        {
         SigmaR[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"SigmaR"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::NO_WRITE
                     ),
                     solidRegions[i],
                     dimensionedScalar("", dimless/dimLength, 0)
             )
         );
        }

    Info<< "Region "<< i<< " Adding fission cross section \n" << endl;
    SigmaF.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(SigmaF[i], j)
        {
         SigmaF[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"SigmaF"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                     solidRegions[i],
                     SigmaF_ref[i][j]
                     
             )
         );
        }

    Info<< "Region "<< i<< " Adding scatter cross section \n" << endl;
    SigmaS.set(i,new PtrList<PtrList<volScalarField>>(energyGroups));
    forAll(SigmaS[i], j)
    {
        SigmaS[i].set(j,new PtrList<volScalarField>(energyGroups));
        forAll(SigmaS[i][j],g)
        {
           SigmaS[i][j].set(
                g,
                new volScalarField(
                    IOobject
                    (
                    	"SigmaS"+Foam::name(j+1)+Foam::name(g+1),
                    	runTime.timeName(),
                    	solidRegions[i],
                    	IOobject::NO_READ,
                    	IOobject::NO_WRITE
                    ),
                solidRegions[i],
                SigmaS_ref[i][j][g]
                )
           );
        }     
    }

    Info<< "Region "<< i << " Adding Qf \n" << endl;
    Qf.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(Qf[i], j)
        {
         Qf[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"Qf"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    //Qf_ref[i][j]
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0.0)
             )
         );
        }

    Info<< "Region "<< i<< " Adding Prescursor density \n" << endl;
    Precursor.set(i,new PtrList<volScalarField>(precGroups)); //
    forAll(Precursor[i], j)
        {
         Precursor[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"Precursor"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::MUST_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i]
             )
         );
        }

    Info<< "Region "<< i<< " Adding prescursor lambda \n" << endl;
    precLambda.set(i,new PtrList<volScalarField>(precGroups)); //
    forAll(precLambda[i], k)
        {
         precLambda[i].set(
             k,
             new volScalarField(
                 IOobject
                     (
                     	"precLambda"+Foam::name(k+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::NO_WRITE
                     ),
                    solidRegions[i],
                    precLambda_ref[i][k]

             )
         );
        }

    Info<< "Region "<< i<< " Adding prescursor beta \n" << endl;
    precBeta.set(i,new PtrList<volScalarField>(precGroups)); //
    forAll(precBeta[i], k)
        {
         precBeta[i].set(
             k,
             new volScalarField(
                 IOobject
                     (
                     	"precBeta"+Foam::name(k+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::NO_WRITE
                     ),
                    solidRegions[i],
                    precBeta_ref[i][k]
             )
         );
        }

    
    
    Info<< "Region "<< i<< " Adding prompt neutron fraction chi promote\n" << endl;
    ChiPrompt.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(ChiPrompt[i], j)
        {
         ChiPrompt[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"ChiPrompt"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                     ),
                    solidRegions[i],
                    chiPrompt_ref[i][j]
             )
         );
        }

    //SigmaS.set(i,new PtrList<PtrList<volScalarField>>(energyGroups));
    //forAll(SigmaS[i], j)
    //{
    //    SigmaS[i].set(j,new PtrList<volScalarField>(energyGroups));
    //    forAll(SigmaS[i][j],g)
    //    {
    //       SigmaS[i][j].set(
    //            g,
    //            new volScalarField(
    //                IOobject
    //                (
    //                	"SigmaS"+Foam::name(j+1)+Foam::name(g+1),
    //                	runTime.timeName(),
    //                	solidRegions[i],
    //                	IOobject::NO_READ,
    //                	IOobject::NO_WRITE
    //                ),
    //            solidRegions[i],
    //            SigmaS_ref[i][j][g]
    //            )
    //       );
    //    }     
    //}
    
    
    Info<< "Region "<< i<< " Adding delayed neutron fraction chi delayed \n" << endl;
    ChiDelayed.set(i,new PtrList<PtrList<volScalarField>>(precGroups));
    forAll(ChiDelayed[i], k)
    {
        ChiDelayed[i].set(k,new PtrList<volScalarField>(energyGroups));
        forAll(ChiDelayed[i][k],g)
        {
           ChiDelayed[i][k].set(
                g,
                new volScalarField(
                    IOobject
                    (
                    	"ChiDelayed"+Foam::name(k+1)+Foam::name(g+1),
                    	runTime.timeName(),
                    	solidRegions[i],
                    	IOobject::NO_READ,
                    	IOobject::NO_WRITE
                    ),
                solidRegions[i],
                chiDelayed_ref[i][k][g]
                )
           );
        }     
    }

    Info<< "Region "<< i<< " Adding fission energy per time Ef \n" << endl;
    Ef.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(Ef[i], j)
        {
         Ef[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"Ef"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::NO_WRITE
                     ),
                    solidRegions[i],
                    Ef_ref[i][j]
             )
         );
        }
    
    Info<< " Adding nu * signmaF\n" << endl;
    NuSigmaF.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(NuSigmaF[i], j)
        {
         NuSigmaF[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"NuSigmaF"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    NuSigmaF_ref[i][j]
             )
         );
        }

    NeutronSource.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(NeutronSource[i], j)
        {
         NeutronSource[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"NeutronSource"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0)
             )
         );
        }

    NeutronSourceAdjoint.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(NeutronSourceAdjoint[i], j)
        {
         NeutronSourceAdjoint[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"NeutronSourceAdjoint"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0)
             )
         );
        }


    GdelayedNeutronSource.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(GdelayedNeutronSource[i], j)
        {
         GdelayedNeutronSource[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"GdelayedNeutronSource"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0)
             )
         );
        }


    ScatteringTerm.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(ScatteringTerm[i], j)
        {
         ScatteringTerm[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"ScatteringTerm"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0)
             )
         );
        }

    ScatteringTermAdjoint.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(ScatteringTermAdjoint[i], j)
        {
         ScatteringTermAdjoint[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"ScatteringTermAdjoint"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0)
             )
         );
        }

    delayedNeutronSource.set(i,new PtrList<volScalarField>(precGroups));
    forAll(delayedNeutronSource[i], j)
        {
         delayedNeutronSource[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"delayedNeutronSource"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::AUTO_WRITE
                     ),
                    solidRegions[i],
                    dimensionedScalar("", dimless/dimVolume/dimTime, 0)
             )
         );
        }


}
