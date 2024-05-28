//1D-Rod Ejection
{

    Info<< "Update Cross Section Started!" << endl;
    //calculate control rod hight
    scalar CRAz=0;
    scalar CRBz=0;
    scalar AxisMeshHight=0.01;

    if (runTime.value()<=4.0 && runTime.value()>0)
    {
        //CRAz=dimensionedScalar("",dimLength,1+runTime.value()*0.03);
        CRAz=1.25+runTime.value()*0.25;
    }
    else if (runTime.value()>4.0)
    {
        //CRAz=dimensionedScalar("",dimLength,1.8);
        CRAz=2.25-(runTime.value()-4.0)*0.25;
    }

    List<scalar> Dg_Mat1(energyGroups,0.0); //without CR
    List<scalar> Dg_Mat2(energyGroups,0.0); //with    CR

    List<scalar> SigmaA_Mat1(energyGroups,0.0);
    List<scalar> SigmaA_Mat2(energyGroups,0.0);

    List<scalar> SigmaF_Mat1(energyGroups,0.0);
    List<scalar> SigmaF_Mat2(energyGroups,0.0);

    List<scalar> NuSigmaF_Mat1(energyGroups,0.0);
    List<scalar> NuSigmaF_Mat2(energyGroups,0.0);

    List<List<scalar>> SigmaS_Mat1(energyGroups,List<scalar>(energyGroups,0.0));
    List<List<scalar>> SigmaS_Mat2(energyGroups,List<scalar>(energyGroups,0.0));

    forAll(solidRegions,i)
    {
        forAll(phi[i],j)
        {
            if (solidRegions[i].name() == "unrodded")
            {
                Dg_Mat1[j]=D_ref[i][j].value();
                SigmaA_Mat1[j]=SigmaA_ref[i][j].value();
                SigmaF_Mat1[j]=SigmaF_ref[i][j].value();
                NuSigmaF_Mat1[j]=NuSigmaF_ref[i][j].value();
                forAll(phi[i],g)
                {
                    SigmaS_Mat1[j][g]=SigmaS_ref[i][j][g].value();
                }
            }

            else if (solidRegions[i].name() == "rodded")
            {
                Dg_Mat2[j]=D_ref[i][j].value();
                SigmaA_Mat2[j]=SigmaA_ref[i][j].value();
                SigmaF_Mat2[j]=SigmaF_ref[i][j].value();
                NuSigmaF_Mat2[j]=NuSigmaF_ref[i][j].value();
                forAll(phi[i],g)
                {
                    SigmaS_Mat2[j][g]=SigmaS_ref[i][j][g].value();
                }
            }

            else
            {
                continue;
            }
            
        }
    }

    //Dg_Mat1={0.01423913 , 0.00356306};
    //Dg_Mat2={0.01423913 , 0.00356306};
    //
    //
    //
    //SigmaT_Mat1={ 1.040206     , 8.766217     };
    //SigmaT_Mat2={ 1.095206     , 9.146217 };
    //
    //
    //
    //SigmaF_Mat1={ 0.6477691      , 11.27328    };
    //SigmaF_Mat2={ 0.6477691      , 11.27328    };
    //
    //
    //
    //SigmaS_Mat1={{ -1.75555 , 0 },
    //             {  1.75555 , 0 }};
    //SigmaS_Mat2={{ -1.75555 , 0 }, 
    //             {  1.75555 , 0 }};
    //
    //

    forAll(solidRegions, i)
    {

        if (solidRegions[i].name()=="unrodded" || solidRegions[i].name()=="rodded") //A rod group lift
        {
            //const auto& NodeCoordList=solidRegions[i].C().internalField();
            const volVectorField& CellCenterField = solidRegions[i].C(); 
            const volScalarField& CellZcoord = CellCenterField.component(2);
            const label& nCells = solidRegions[i].nCells(); 
            Info<<"The total cell numbers: "<<nCells<<endl;
            
            forAll(phi[i], j)
            {
                //const auto& fieldNodes = solidRegions[i].Cf();

                label nPhis = 0;
                
                forAll(phi[i][j], nodeIndex)
                {
                    nPhis += 1;
                    if ((-CellZcoord[nodeIndex]+CRAz) < AxisMeshHight/2 && (-CellZcoord[nodeIndex]+CRAz) > -AxisMeshHight/2)
                    {
                        scalar hCR=CRAz-CellZcoord[nodeIndex] + AxisMeshHight/2;
                        SigmaF[i][j][nodeIndex]=hCR / AxisMeshHight * SigmaF_Mat1[j] + (AxisMeshHight - hCR )/ AxisMeshHight * SigmaF_Mat2[j];
                        NuSigmaF[i][j][nodeIndex]=hCR / AxisMeshHight * NuSigmaF_Mat1[j] + (AxisMeshHight - hCR )/ AxisMeshHight * NuSigmaF_Mat2[j];
                    }
                    else if ((-CellZcoord[nodeIndex]+CRAz) <= -AxisMeshHight/2) //with-CR
                    {
                        SigmaF[i][j][nodeIndex]= SigmaF_Mat2[j];
                        NuSigmaF[i][j][nodeIndex]= NuSigmaF_Mat2[j];
                    }
                    else if ((-CellZcoord[nodeIndex]+CRAz) >= AxisMeshHight/2) //without-CR
                    {
                        SigmaF[i][j][nodeIndex]= SigmaF_Mat1[j];
                        NuSigmaF[i][j][nodeIndex]= NuSigmaF_Mat1[j];
                    }

                }
                Info<<"The total phi numbers: "<<nPhis<<endl;
            }
        }

        else
        {
            continue;
        }            
    }
    
    Info<< "Update Cross Section End!" << endl;
}