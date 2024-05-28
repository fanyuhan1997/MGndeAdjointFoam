IOdictionary nuclearPropertiesDictInit	//Define the reference dictionary to read the number of neutron groups, precursor groups and decay pr
(
    IOobject
    (
        "nuclearProperties",
        runTime.constant(),
        solidRegions[0],
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);

const label energyGroups(readInt(nuclearPropertiesDictInit.lookup("energyGroups")));
const label precGroups(readInt(nuclearPropertiesDictInit.lookup("precGroups")));
const scalar KeffIni(readScalar(nuclearPropertiesDictInit.lookup("KeffIni")));
const scalar Ktol(readScalar(nuclearPropertiesDictInit.lookup("Ktol"))); 
const scalar Powertol(readScalar(nuclearPropertiesDictInit.lookup("Powertol"))); 

const label nRegions = solidRegions.size();
//
PtrList<PtrList<dimensionedScalar>> D_ref(nRegions);
//PtrList<PtrList<dimensionedScalar>> SigmaT_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> SigmaA_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> NuSigmaF_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> SigmaF_ref(nRegions);
PtrList<PtrList<PtrList<dimensionedScalar>>> SigmaS_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> InverseVelocity_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> Nu_ref(nRegions);

PtrList<PtrList<dimensionedScalar>> Ef_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> chiPrompt_ref(nRegions);
PtrList<PtrList<PtrList<dimensionedScalar>>> chiDelayed_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> Qf_ref(nRegions);
PtrList<dimensionedScalar> Beta_ref(nRegions);

//PtrList<PtrList<dimensionedScalar>> Precursor_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> precLambda_ref(nRegions);
PtrList<PtrList<dimensionedScalar>> precBeta_ref(nRegions);

// Now we need nprecGroups*nRegions for what regards precEqns, despite of precBetaTot that is a number per each region.

//PtrList<PtrList<dimensionedScalar>> precLambda_ref(nRegions);
//PtrList<PtrList<dimensionedScalar>> precBeta_ref(nRegions);
//PtrList<dimensionedScalar> precBetaTot_ref(nRegions);

// Now we need nDecGroupsInit*nRegions for what regards decEqns, despite of decBetaTot that is a single number per each region.

//PtrList<PtrList<dimensionedScalar>> decLambda(nRegions);
//PtrList<PtrList<dimensionedScalar>> decBeta(nRegions);
//PtrList<dimensionedScalar> decBetaTot(nRegions);




forAll(solidRegions,i) //j-regions
{

    word iName="_"+name(i + 1);
    IOdictionary nuclearPropertiesDict	//Define the dictionary
    (
		IOobject
      (
        "nuclearProperties",
        runTime.constant(),
        solidRegions[i],
	      IOobject::MUST_READ_IF_MODIFIED,
	      IOobject::NO_WRITE
      )
    );

    dictionary neutronTransportDict(nuclearPropertiesDict.subDict("neutronTransport"));

    List<scalar> D_ref_(neutronTransportDict.lookup("DRef"));//

    //List<scalar> SigmaT_ref_(neutronTransportDict.lookup("SigmaTRef"));/

    List<scalar> SigmaA_ref_(neutronTransportDict.lookup("SigmaARef"));//

    List<scalar> NuSigmaF_ref_(neutronTransportDict.lookup("NuSigmaFRef")); //
    List<scalar> SigmaF_ref_(neutronTransportDict.lookup("SigmaFRef")); //

    List<List<scalar>> SigmaS_ref_(neutronTransportDict.lookup("SigmaSRef")); //

    List<scalar> InverseVelocity_(neutronTransportDict.lookup("InverseVelocity"));//

    List<scalar> Nu_(neutronTransportDict.lookup("Nu"));//

    

    List<scalar> Ef_(neutronTransportDict.lookup("Ef"));//

    List<scalar> chiPrompt_(neutronTransportDict.lookup("chiPrompt"));//

    List<List<scalar>> chiDelayed_(neutronTransportDict.lookup("chiDelayed"));//

    const scalar Beta_(readScalar(neutronTransportDict.lookup("Beta"))); //



    //List<scalar> Precursor_(neutronTransportDict.lookup("Precursor"));/

    List<scalar> precLambda_(neutronTransportDict.lookup("precLambda"));//

    List<scalar> precBeta_(neutronTransportDict.lookup("precBeta"));//



    //List<scalar> Qf_(neutronTransportDict.lookup("Qf"));/

    


    D_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    //SigmaT_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    SigmaA_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    NuSigmaF_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    SigmaF_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    SigmaS_ref.set(i,new PtrList<PtrList<dimensionedScalar>> (energyGroups));
    InverseVelocity_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    Nu_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    
    Ef_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));
    chiPrompt_ref.set(i,new PtrList<dimensionedScalar> (energyGroups));

    chiDelayed_ref.set(i,new PtrList<PtrList<dimensionedScalar>> (precGroups));

    Beta_ref.set(i,new dimensionedScalar("Beta_ref"+iName,dimless,Beta_));

    //Precursor_ref.set(i,new PtrList<dimensionedScalar>(precGroups));
    precLambda_ref.set(i,new PtrList<dimensionedScalar>(precGroups));
    precBeta_ref.set(i,new PtrList<dimensionedScalar>(precGroups));

    //Beta_ref.set(i,new PtrList<dimensionedScalar>(energyGroups));
    //Qf_ref.set(i,new PtrList<dimensionedScalar>(energyGroups));

    
    forAll(precLambda_,k)
    {

      word ikName="_" + name(i + 1)  + name(k + 1);
      //Precursor_ref[i].set(k, new dimensionedScalar("Precursor_ref"   + ikName , dimless/dimVolume/dimTime, Precursor_[k]));
      precLambda_ref[i].set(k, new dimensionedScalar("precLambda_ref" + ikName , dimless/dimTime, precLambda_[k]));
      precBeta_ref[i].set(k, new dimensionedScalar("precBeta_ref"   + ikName , dimless, precBeta_[k]));

      chiDelayed_ref[i].set(k, new PtrList<dimensionedScalar> (energyGroups));
      forAll(D_ref_,g) //
      {
        word ikgName = "_"+name(i + 1)  + name(k + 1) + name(g + 1);
        chiDelayed_ref[i][k].set(g, new dimensionedScalar("chiDelayed_ref" + ikgName , dimless, chiDelayed_[k][g]));
        Info<<"region "<< i <<" Prcursor "<< k << " Group "<< g <<" chiDelayed_ref = "<<chiDelayed_ref[i][k][g]<<endl;
      }


    }
    
    forAll(D_ref_, j) //j
    {
      word ijName = "_"+name(i + 1)  + name(j + 1);

      D_ref[i].set(j, new dimensionedScalar("D_ref" + ijName , dimless*dimLength, D_ref_[j]));
      //SigmaT_ref[i].set(j, new dimensionedScalar("SigmaT_ref" + ijName , dimless/dimLength, SigmaT_ref_[j]));
      SigmaA_ref[i].set(j, new dimensionedScalar("SigmaA_ref" + ijName , dimless/dimLength, SigmaA_ref_[j]));
      NuSigmaF_ref[i].set(j, new dimensionedScalar("NuSigmaF_ref" + ijName , dimless/dimLength, NuSigmaF_ref_[j]));
      SigmaF_ref[i].set(j, new dimensionedScalar("SigmaF_ref" + ijName , dimless/dimLength, SigmaF_ref_[j]));


      SigmaS_ref[i].set(j,new PtrList<dimensionedScalar>(energyGroups));
      forAll(D_ref_,k) //
      {
        word ijkName = "_"+name(i + 1)  + name(j + 1) + name(k + 1);
        SigmaS_ref[i][j].set(k, new dimensionedScalar("SigmaS_ref" + ijkName , dimless/dimLength, SigmaS_ref_[j][k]));
        Info<<"region "<<i<<" group "<<j<<"SigmaS_ref = "<<SigmaS_ref[i][j][k]<<endl;
      }

      Nu_ref[i].set(j, new dimensionedScalar("Nu_ref" + ijName , dimless, Nu_[j]));
      
      Ef_ref[i].set(j, new dimensionedScalar("Ef_ref" + ijName , dimless, Ef_[j]));
      InverseVelocity_ref[i].set(j, new dimensionedScalar("InverseVelocity_ref" + ijName , dimless/dimLength*dimTime, InverseVelocity_[j]));
      chiPrompt_ref[i].set(j, new dimensionedScalar("chiPrompt_ref" + ijName , dimless , chiPrompt_[j]));

      //chiDelayed_ref[i].set(j, new dimensionedScalar("chiDelayed_ref" + ijName , dimless , chiDelayed_[j]));



      
      Info<<"region "<<i<<" group "<<j<<" NuSigmaF_ref = "<<NuSigmaF_ref[i][j]<<endl;
      
      
      //Qf_ref[i].set(j, new dimensionedScalar("Qf_ref" + ijName , dimless , Qf_[j]));

    }

    

    //dictionary neutronPrecursorsDict(nuclearPropertiesDict.subDict("delayedNeutronPrecursors"));
    //List<scalar> precLambda_ref_(neutronPrecursorsDict.lookup("Lambda"));
    //List<scalar> precBeta_ref_(neutronPrecursorsDict.lookup("Beta"));
    //scalar precBetaTot_ref_ = 0.0;

    //forAll(precLambda_ref_, j)
    //{
    //  word ijName = "_"+name(i + 1)  + name(j + 1);
    //
    //  precLambda_ref[i].set(j, new dimensionedScalar("precLambda" + ijName , dimless/dimTime, precLambda_ref_[j]));
    //  precBeta_ref[i].set(j, new dimensionedScalar("precBeta" + ijName , dimless, precBeta_ref_[j]));
    //  precBetaTot_ref_ += precBeta_ref[i][j].value();
    //
    //}
    //
    //precBetaTot_ref.set(i, new dimensionedScalar("precBetaTot_" + name(i+1), dimless, precBetaTot_ref_ ));
}
Info<< "\n End of readSolidNuclearProperties" << endl;







