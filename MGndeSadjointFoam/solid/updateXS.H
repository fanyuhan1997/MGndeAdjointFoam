//"step reactivity"
{
    forAll(solidRegions, i)
    {
        if (solidRegions[i].name()=="region1")
        {
            forAll(phi[i], j)
            {
                if (j==1 && runTime.value() > 0.0)
                    SigmaT[i][j]=dimensionedScalar("", dimless/dimLength,0.15e2*0.97666);
                else
                {
                    continue;
                }
            }
        }
        else
        {
            continue;
        }            
    }
} 


//"linear reactivity"
//{
//    forAll(solidRegions, i)
//    {
//        if (solidRegions[i].name()=="region1")
//        {
//            forAll(phi[i], j)
//            {
//                if (j==1 && runTime.value()<=0.2)
//                    SigmaT[i][j]=dimensionedScalar("", dimless/dimLength,0.15e2*(1-runTime.value()*0.11667));
//                else if (j==1 && runTime.value()>0.2)
//                    SigmaT[i][j]=dimensionedScalar("", dimless/dimLength,0.15e2*0.97666);
//                    
//                else
//                {
//                    continue;
//                }
//            }
//        }
//        else
//        {
//            continue;
//        }            
//    }
//    
//}