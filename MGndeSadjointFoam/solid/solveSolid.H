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
 


        forAll(phi[i], j) //
          {

          ScatteringTerm[i][j] *= 0;

          forAll(phi[i], g) //updat source term
              {
                ScatteringTerm[i][j]  += SigmaS[i][g][j] * phi[i][g];// 
                       
              }
          
          NeutronSource[i][j]*=0;
          NeutronSource[i][j]= ScatteringTerm[i][j] + (1-Beta[i]) * TotGroupFissions0[i] * ChiPrompt[i][j]  / Keff ; //DelayedTerm; 
          
          GdelayedNeutronSource[i][j] *= 0;
          forAll(Precursor[i], k)
          {
              GdelayedNeutronSource[i][j] += precLambda[i][k] * Precursor[i][k] * ChiDelayed[i][k][j];
          }
          //prevFissions+=TotGroupFissions[i];

          // solve steady diffusion equation
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
        
      
          
    } 

//thermo.correct();

//Info<< "Min/max T:" << min(thermo.T()).value() << ' '
//    << max(thermo.T()).value() << endl;
}