{
    while (pimple.correctNonOrthogonal())
    {

        //Beta[i]*=0;
        //forAll(Precursor[i], k)
        //{
        //  Beta[i]+=precBeta[i][k];
        //}

        forAll(phiAdjoint[i], j) //
          {
          //volScalarField& phiAdjointIJ = phiAdjoint[i][j]; //
          //volScalarField& DgIJ= Dg[i][j];
          //volScalarField& VelIJ= Vel[i][j];
          //volScalarField& SigmaTIJ = SigmaT[i][j]; //

          //volScalarField ScatteringTerm[i][j] = zeroField[i]*dimensionedScalar("", dimless/dimVolume/dimTime, 0);  
          ScatteringTermAdjoint[i][j] *= 0;

          forAll(phiAdjoint[i], g) // update source term
              {
                ScatteringTermAdjoint[i][j]  += SigmaS[i][j][g] * phiAdjoint[i][g];//        
              }
          
          NeutronSourceAdjoint[i][j]*=0;
          NeutronSourceAdjoint[i][j]= ScatteringTermAdjoint[i][j] +  NuSigmaF[i][j] * ChiAdjointFlux[i] / KeffAdjoint ; //DelayedTerm; (1-Beta[i]) * TotGroupFissions0[i] * ChiPrompt[i][j]
          
          //GdelayedNeutronSource[i][j] *= 0;
          //forAll(Precursor[i], k)
          //{
          //    GdelayedNeutronSource[i][j] += precLambda[i][k] * Precursor[i][k] * ChiDelayed[i][k][j];
          //}
          

          // solve adjoint diffusion equation
          fvScalarMatrix phiAdjointEqn
          (
                
              - fvm::laplacian(Dg[i][j], phiAdjoint[i][j])
              + fvm::Sp(SigmaR[i][j],phiAdjoint[i][j])
              == 
                NeutronSourceAdjoint[i][j]
          );

          phiAdjointEqn.relax(); // relax
          //fvConstraints.constrain(phiEqn);
          phiAdjointEqn.solve(); //solve
          //fvConstraints.constrain(phiIJ);        

          } 
        //
      
          
    } 

//thermo.correct();

//Info<< "Min/max T:" << min(thermo.T()).value() << ' '
//    << max(thermo.T()).value() << endl;
}