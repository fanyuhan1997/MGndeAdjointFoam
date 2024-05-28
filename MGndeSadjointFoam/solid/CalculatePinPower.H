std::fstream outfile_PP;
outfile_PP.open("PinPower.dat",std::fstream::out);

scalar PinPower=0;
forAll(solidRegions, i)
{
    PinPower*=0;
    PinPower = gSum(fvc::volumeIntegrate(TotGroupFissionEnergy[i])); 
    Info<< "region: "<< solidRegions[i].name() << " Power="<<PinPower<<endl; 
    if (outfile_PP.is_open())
    {
        outfile_PP<<"region: "<< solidRegions[i].name() << " Power="<<PinPower<<std::endl; 
        outfile_PP.close();
    }  
    else
    {
        outfile_PP.open("PinPower.dat",std::fstream::app);
        outfile_PP<<"region: "<< solidRegions[i].name() << " Power="<<PinPower<<std::endl;
        outfile_PP.close();
    }
}