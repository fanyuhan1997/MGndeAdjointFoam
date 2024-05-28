// Initialise solid field pointer lists
#include "nuclearProperties.H" //
//PtrList<solidThermo> thermoSolid(solidRegions.size()); /
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size()); //
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size()); //

PtrList<PtrList<volScalarField>> phi(solidRegions.size()); //
PtrList<PtrList<volScalarField>> phiPre(solidRegions.size()); //)
PtrList<PtrList<volScalarField>> phiErr(solidRegions.size()); //

PtrList<PtrList<volScalarField>> Dg(solidRegions.size()); //
PtrList<PtrList<volScalarField>> InverseVel(solidRegions.size()); //
PtrList<PtrList<volScalarField>> SigmaA(solidRegions.size()); //
PtrList<PtrList<volScalarField>> SigmaR(solidRegions.size()); //
PtrList<PtrList<volScalarField>> SigmaF(solidRegions.size()); //
PtrList<PtrList<PtrList<volScalarField>>> SigmaS(solidRegions.size());//，g* 
PtrList<PtrList<volScalarField>> Qs(solidRegions.size()); 
PtrList<PtrList<volScalarField>> Qf(solidRegions.size()); 

PtrList<PtrList<volScalarField>> ChiPrompt(solidRegions.size()); //
PtrList<PtrList<PtrList<volScalarField>>> ChiDelayed(solidRegions.size()); //
PtrList<PtrList<volScalarField>> NuSigmaF(solidRegions.size());  //
PtrList<PtrList<volScalarField>> NeutronSource(solidRegions.size()); //
PtrList<PtrList<volScalarField>> ScatteringTerm(solidRegions.size()); //

PtrList<PtrList<volScalarField>> delayedNeutronSource(solidRegions.size()); //
PtrList<PtrList<volScalarField>> GdelayedNeutronSource(solidRegions.size()); //
PtrList<PtrList<volScalarField>> Precursor(solidRegions.size()); //
PtrList<PtrList<volScalarField>> precLambda(solidRegions.size());//
PtrList<PtrList<volScalarField>> precBeta(solidRegions.size());//

PtrList<volScalarField> zeroField(solidRegions.size()); //

PtrList<volScalarField> Beta(solidRegions.size()); //
PtrList<PtrList<volScalarField>> Nu(solidRegions.size()); //

PtrList<volScalarField> NuRfTot(solidRegions.size()); //
PtrList<PtrList<volScalarField>> Ef(solidRegions.size()); //

PtrList<volScalarField> TotGroupFissions(solidRegions.size()); //
PtrList<volScalarField> TotGroupFissions0(solidRegions.size()); //

PtrList<volScalarField> TotGroupFissionEnergy(solidRegions.size()); //

PtrList<PtrList<volScalarField>> phiAdjoint(solidRegions.size()); 
PtrList<PtrList<volScalarField>> lapphi(solidRegions.size()); 
PtrList<PtrList<volScalarField>> DgZero(solidRegions.size()); //
PtrList<PtrList<volScalarField>> SigmaAZero(solidRegions.size()); //
PtrList<PtrList<PtrList<volScalarField>>> SigmaSZero(solidRegions.size());//g* 
PtrList<PtrList<volScalarField>> NuSigmaFZero(solidRegions.size());  //

PtrList<volScalarField> RhoPerturbation0(solidRegions.size()); //Numerator of reactivity
PtrList<volScalarField> RhoPerturbation1(solidRegions.size()); //Denominator of reactivity

dimensionedScalar Keff("Keff", dimless, scalar(KeffIni));
dimensionedScalar RhoNormal("RhoNormal", dimless, scalar(0.0));
dimensionedScalar MiNiNum("", dimless/dimArea/dimTime, 1e-16);
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

scalar TotalFissionPower=1.0;
scalar TotalFissionPower0=1.0;
scalar TotalFissionPower1=1.0;

scalar PerturbationFissions0=1.0;
scalar PerturbationFissions1=1.0;

scalar Kerr=KeffIni;
scalar MaxPhiErr=1.0;
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
        dimensionedScalar("",dimless, 0.0)
        //Beta_ref[i]
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
        forAll(phi[i], j) //j
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

    Info<< "Region "<< i<< " Adding flux Previous \n" << endl;
    phiPre.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(phiPre[i], j) //j
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
        forAll(phiErr[i], j) //j
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

    DgZero.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(DgZero[i], j)
        {
            DgZero[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"DgZero"+Foam::name(j+1),
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
                     	IOobject::NO_WRITE
                     ),
                     solidRegions[i],
                     SigmaA_ref[i][j]
             )
         );
        }

    SigmaAZero.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(SigmaAZero[i], j)
        {
         SigmaAZero[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"SigmaAZero"+Foam::name(j+1),
                     	runTime.timeName(),
                     	solidRegions[i],
                     	IOobject::NO_READ,
                     	IOobject::NO_WRITE
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

    Info<< "Region "<< i<< " Adding scatter cross section \n" << endl;
    SigmaSZero.set(i,new PtrList<PtrList<volScalarField>>(energyGroups));
    forAll(SigmaSZero[i], j)
    {
        SigmaSZero[i].set(j,new PtrList<volScalarField>(energyGroups));
        forAll(SigmaSZero[i][j],g)
        {
           SigmaSZero[i][j].set(
                g,
                new volScalarField(
                    IOobject
                    (
                    	"SigmaSZero"+Foam::name(j+1)+Foam::name(g+1),
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
    
    Info<< " Adding nu * signmaF\n" << endl;
    NuSigmaFZero.set(i,new PtrList<volScalarField>(energyGroups));
    forAll(NuSigmaFZero[i], j)
        {
         NuSigmaFZero[i].set(
             j,
             new volScalarField(
                 IOobject
                     (
                     	"NuSigmaFZero"+Foam::name(j+1),
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

    Info<< "Region "<< i<< " Adding flux \n" << endl;
    phiAdjoint.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(phiAdjoint[i], j) //j
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

    Info<< "Region "<< i<< " Adding flux \n" << endl;
    lapphi.set(i,new PtrList<volScalarField>(energyGroups));
        forAll(lapphi[i], j) //j-
        {
            lapphi[i].set(
                j,
                new volScalarField(
                    IOobject
                        (
                        	"lapphi"+Foam::name(j+1),
                        	runTime.timeName(),
                        	solidRegions[i],
                        	IOobject::NO_READ, //
                        	IOobject::NO_WRITE 
                        ),
                        solidRegions[i],
                        dimensionedScalar("", dimless/dimArea/dimArea/dimTime,0)
                )
            );
        }


    Info<< " ******" << endl;
    RhoPerturbation0.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"RhoPerturbation0",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::NO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("", dimless/dimLength/dimArea/dimArea/dimTime/dimTime, 0)
    	)
    );

    Info<< " ******" << endl;
    RhoPerturbation1.set
    (
    	i,
    	new volScalarField
    	(
    		IOobject
    		(
    			"RhoPerturbation1",
    			runTime.timeName(),
    			solidRegions[i],
    			IOobject::NO_READ,
    			IOobject::NO_WRITE
    		),
    	solidRegions[i],
    	dimensionedScalar("", dimless/dimLength/dimArea/dimArea/dimTime/dimTime, 0)
    	)
    );
}
