{

    
    while (pimple.correctNonOrthogonal())
    {
        //fvScalarMatrix eEqn
        //(
        //    fvm::ddt(rho, e)
        //  + thermo.divq(e)
        //  ==
        //    fvModels.source(rho, e)
        //);
        //eEqn.relax(); // relax
        //fvConstraints.constrain(eEqn);
        //eEqn.solve(); //solve
        //fvConstraints.constrain(e);

        //neutron  
        //


        volScalarField FissionTerm = zeroField[i]*dimensionedScalar("", dimless/dimVolume/dimTime, 0);
        forAll(phi[i], g) 
        {
          FissionTerm +=  (1-Beta[i]) * NuSigmaF[i][g] *phi[i][g] ;    
        }

        Beta[i]*=0;
        forAll(Precursor[i], k)
        {
          Beta[i]+=precBeta[i][k];
        }

        forAll(phi[i], j) //i,j
          {
          //volScalarField& phiIJ = phi[i][j]; /phiIJ
          //volScalarField& DgIJ= Dg[i][j];
          //volScalarField& VelIJ= Vel[i][j];
          //volScalarField& SigmaTIJ = SigmaT[i][j]; /

          
          //volScalarField ScatteringTerm[i][j] = zeroField[i]*dimensionedScalar("", dimless/dimVolume/dimTime, 0);  
          ScatteringTerm[i][j] *= 0;

          forAll(phi[i], g) 
              {
                ScatteringTerm[i][j]  += SigmaS[i][g][j] * phi[i][g];//        
                Qf[i][g]        =   NuSigmaF[i][g] *phi[i][g];
              }
          
          NeutronSource[i][j]*=0;
          NeutronSource[i][j]= ScatteringTerm[i][j] + (1-Beta[i]) * TotGroupFissions0[i] * ChiPrompt[i][j]  / Keff ; //DelayedTerm; 
          
          GdelayedNeutronSource[i][j] *= 0;
          forAll(Precursor[i], k)
          {
              GdelayedNeutronSource[i][j] += precLambda[i][k] * Precursor[i][k] * ChiDelayed[i][k][j];
          }
          //prevFissions+=TotGroupFissions[i];

          //
          fvScalarMatrix phiEqn
          (
                fvm::ddt(InverseVel[i][j],phi[i][j])
              - fvm::laplacian(Dg[i][j], phi[i][j])
              + fvm::Sp(SigmaR[i][j],phi[i][j])
              == 
                NeutronSource[i][j]
              + GdelayedNeutronSource[i][j]
          );

          phiEqn.relax(); // relax
          //fvConstraints.constrain(phiEqn);
          phiEqn.solve(); //solve
          //fvConstraints.constrain(phiIJ);        

          } 
        
      
          
    } //pimpl

//thermo.correct();/correct(

//Info<< "Min/max T:" << min(thermo.T()).value() << ' '
//    << max(thermo.T()).value() << endl;
}