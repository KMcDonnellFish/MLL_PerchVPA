#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <Perch2025.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  fyear.allocate("fyear");
  lyearread.allocate("lyearread");
  lyear.allocate("lyear");
  nages.allocate("nages");
  ReadInObsKill.allocate(fyear,lyearread,1,nages,"ReadInObsKill");
  ReadInObsGillNet.allocate(fyear,lyearread+1,1,nages,"ReadInObsGillNet");
  ReadInObsOffshoreGillNet.allocate(1999,lyearread+1,1,nages,"ReadInObsOffshoreGillNet");
  ReadInObsTrawlYOY.allocate(fyear,lyearread+2,"ReadInObsTrawlYOY");
  ReadInObsTrawlOnePlus.allocate(1989,lyearread+1,"ReadInObsTrawlOnePlus");
  M.allocate(1,nages,"M");
  ReadInObsEffort.allocate(fyear,lyearread,"ReadInObsEffort");
  KillSigma.allocate("KillSigma");
  GillNetSigma.allocate("GillNetSigma");
  OffshoreGillNetSigma.allocate("OffshoreGillNetSigma");
  TrawlYOYSigma.allocate("TrawlYOYSigma");
  TrawlOnePlusSigma.allocate("TrawlOnePlusSigma");
  FSigma.allocate("FSigma");
  KillSampSize.allocate("KillSampSize");
  GillNetSampSize.allocate("GillNetSampSize");
  OffshoreGillNetSampSize.allocate("OffshoreGillNetSampSize");
}

void model_parameters::initializationfunction(void)
{
  LogNAgeOne.set_initial_value(14.5);
  FByYear.set_initial_value(0.5);
  EstKillGammaSelectivityShape.set_initial_value(4.0);
  EstKillGammaSelectivityScale.set_initial_value(0.7);
  EstGillNetGammaSelectivityShape.set_initial_value(2.0);
  EstGillNetGammaSelectivityScale.set_initial_value(0.2);
  EstOffshoreGillNetGammaSelectivityShape.set_initial_value(2.0);
  EstOffshoreGillNetGammaSelectivityScale.set_initial_value(0.2);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  LogNAgeOne.allocate(fyear-nages+1,lyear+1,7.6,20,1,"LogNAgeOne");
  FByYear.allocate(fyear,lyear,0.01,2,2,"FByYear");
  LogEffortq.allocate(-25,-2,1,"LogEffortq");
  EstKillSelectivity14.allocate(1,4,0,5,-3,"EstKillSelectivity14");
  EstKillSelectivity6.allocate(0,5,-3,"EstKillSelectivity6");
  EstKillGammaSelectivityShape.allocate(1.1,10,3,"EstKillGammaSelectivityShape");
  EstKillGammaSelectivityScale.allocate(0.001,5,3,"EstKillGammaSelectivityScale");
  EstKillGammaSelectivityLocation.allocate(-1.0,1.0,-3,"EstKillGammaSelectivityLocation");
  LogGillNetq.allocate(-25,-2,1,"LogGillNetq");
  EstGillNetSelectivity13.allocate(1,3,0,5,-3,"EstGillNetSelectivity13");
  EstGillNetSelectivity56.allocate(5,6,0,5,-3,"EstGillNetSelectivity56");
  EstGillNetSelectivityParameter1.allocate(-5,10,-3,"EstGillNetSelectivityParameter1");
  EstGillNetSelectivityParameter2.allocate(-2.5,2.5,-3,"EstGillNetSelectivityParameter2");
  EstGillNetSelectivityParameter3.allocate(-5,10,-3,"EstGillNetSelectivityParameter3");
  EstGillNetSelectivityParameter4.allocate(-2.5,2.5,-3,"EstGillNetSelectivityParameter4");
  EstGillNetGammaSelectivityShape.allocate(1.1,10,3,"EstGillNetGammaSelectivityShape");
  EstGillNetGammaSelectivityScale.allocate(0.001,5,3,"EstGillNetGammaSelectivityScale");
  EstGillNetGammaSelectivityLocation.allocate(-1.0,1.0,-3,"EstGillNetGammaSelectivityLocation");
  LogOffshoreGillNetq.allocate(-25,-2,1,"LogOffshoreGillNetq");
  EstOffshoreGillNetSelectivity13.allocate(1,3,0,5,-3,"EstOffshoreGillNetSelectivity13");
  EstOffshoreGillNetSelectivity56.allocate(5,6,0,5,-3,"EstOffshoreGillNetSelectivity56");
  EstOffshoreGillNetSelectivityParameter1.allocate(-5,10,-3,"EstOffshoreGillNetSelectivityParameter1");
  EstOffshoreGillNetSelectivityParameter2.allocate(-2.5,2.5,-3,"EstOffshoreGillNetSelectivityParameter2");
  EstOffshoreGillNetSelectivityParameter3.allocate(-5,10,-3,"EstOffshoreGillNetSelectivityParameter3");
  EstOffshoreGillNetSelectivityParameter4.allocate(-2.5,2.5,-3,"EstOffshoreGillNetSelectivityParameter4");
  EstOffshoreGillNetGammaSelectivityShape.allocate(1.1,10,3,"EstOffshoreGillNetGammaSelectivityShape");
  EstOffshoreGillNetGammaSelectivityScale.allocate(0.001,5,3,"EstOffshoreGillNetGammaSelectivityScale");
  EstOffshoreGillNetGammaSelectivityLocation.allocate(-1.0,1.0,-3,"EstOffshoreGillNetGammaSelectivityLocation");
  LogTrawlYOYq.allocate(-25,-2,1,"LogTrawlYOYq");
  LogTrawlOnePlusq.allocate(-25,-2,1,"LogTrawlOnePlusq");
  NAgeOne.allocate(fyear-nages+1,lyear+1,"NAgeOne");
  #ifndef NO_AD_INITIALIZE
    NAgeOne.initialize();
  #endif
  Effortq.allocate("Effortq");
  #ifndef NO_AD_INITIALIZE
  Effortq.initialize();
  #endif
  GillNetq.allocate("GillNetq");
  #ifndef NO_AD_INITIALIZE
  GillNetq.initialize();
  #endif
  OffshoreGillNetq.allocate("OffshoreGillNetq");
  #ifndef NO_AD_INITIALIZE
  OffshoreGillNetq.initialize();
  #endif
  TrawlYOYq.allocate("TrawlYOYq");
  #ifndef NO_AD_INITIALIZE
  TrawlYOYq.initialize();
  #endif
  TrawlOnePlusq.allocate("TrawlOnePlusq");
  #ifndef NO_AD_INITIALIZE
  TrawlOnePlusq.initialize();
  #endif
  ObsKill.allocate(fyear,lyear,1,nages,"ObsKill");
  #ifndef NO_AD_INITIALIZE
    ObsKill.initialize();
  #endif
  ObsGillNet.allocate(fyear,lyear+1,1,nages,"ObsGillNet");
  #ifndef NO_AD_INITIALIZE
    ObsGillNet.initialize();
  #endif
  ObsOffshoreGillNet.allocate(1999,lyear+1,1,nages,"ObsOffshoreGillNet");
  #ifndef NO_AD_INITIALIZE
    ObsOffshoreGillNet.initialize();
  #endif
  ObsTrawlYOY.allocate(fyear,lyear+1,"ObsTrawlYOY");
  #ifndef NO_AD_INITIALIZE
    ObsTrawlYOY.initialize();
  #endif
  ObsTrawlOnePlus.allocate(1989,lyear+1,"ObsTrawlOnePlus");
  #ifndef NO_AD_INITIALIZE
    ObsTrawlOnePlus.initialize();
  #endif
  ObsEffort.allocate(fyear,lyear,"ObsEffort");
  #ifndef NO_AD_INITIALIZE
    ObsEffort.initialize();
  #endif
  ObsKillByYear.allocate(fyear,lyear,"ObsKillByYear");
  #ifndef NO_AD_INITIALIZE
    ObsKillByYear.initialize();
  #endif
  ObsGillNetByYear.allocate(fyear,lyear+1,"ObsGillNetByYear");
  #ifndef NO_AD_INITIALIZE
    ObsGillNetByYear.initialize();
  #endif
  ObsOffshoreGillNetByYear.allocate(1999,lyear+1,"ObsOffshoreGillNetByYear");
  #ifndef NO_AD_INITIALIZE
    ObsOffshoreGillNetByYear.initialize();
  #endif
  ObsKillPropAge.allocate(fyear,lyear,1,nages,"ObsKillPropAge");
  #ifndef NO_AD_INITIALIZE
    ObsKillPropAge.initialize();
  #endif
  ObsGillNetPropAge.allocate(fyear,lyear+1,1,nages,"ObsGillNetPropAge");
  #ifndef NO_AD_INITIALIZE
    ObsGillNetPropAge.initialize();
  #endif
  ObsOffshoreGillNetPropAge.allocate(1999,lyear+1,1,nages,"ObsOffshoreGillNetPropAge");
  #ifndef NO_AD_INITIALIZE
    ObsOffshoreGillNetPropAge.initialize();
  #endif
  PredKillByYear.allocate(fyear,lyear,"PredKillByYear");
  #ifndef NO_AD_INITIALIZE
    PredKillByYear.initialize();
  #endif
  PredGillNetByYear.allocate(fyear,lyear+1,"PredGillNetByYear");
  #ifndef NO_AD_INITIALIZE
    PredGillNetByYear.initialize();
  #endif
  PredOffshoreGillNetByYear.allocate(1999,lyear+1,"PredOffshoreGillNetByYear");
  #ifndef NO_AD_INITIALIZE
    PredOffshoreGillNetByYear.initialize();
  #endif
  PredKillPropAge.allocate(fyear,lyear,1,nages,"PredKillPropAge");
  #ifndef NO_AD_INITIALIZE
    PredKillPropAge.initialize();
  #endif
  PredGillNetPropAge.allocate(fyear,lyear+1,1,nages,"PredGillNetPropAge");
  #ifndef NO_AD_INITIALIZE
    PredGillNetPropAge.initialize();
  #endif
  PredOffshoreGillNetPropAge.allocate(1999,lyear+1,1,nages,"PredOffshoreGillNetPropAge");
  #ifndef NO_AD_INITIALIZE
    PredOffshoreGillNetPropAge.initialize();
  #endif
  OnePlusAbundanceByYear.allocate(fyear,lyear+1,"OnePlusAbundanceByYear");
  #ifndef NO_AD_INITIALIZE
    OnePlusAbundanceByYear.initialize();
  #endif
  ThreePlusAbundanceByYear.allocate(fyear,lyear+1,"ThreePlusAbundanceByYear");
  #ifndef NO_AD_INITIALIZE
    ThreePlusAbundanceByYear.initialize();
  #endif
  MaxAgeMinusLocation.allocate(1,nages,"MaxAgeMinusLocation");
  #ifndef NO_AD_INITIALIZE
    MaxAgeMinusLocation.initialize();
  #endif
  InstantaneousF.allocate(fyear,lyear,1,nages,"InstantaneousF");
  #ifndef NO_AD_INITIALIZE
    InstantaneousF.initialize();
  #endif
  InstantaneousZ.allocate(fyear,lyear,1,nages,"InstantaneousZ");
  #ifndef NO_AD_INITIALIZE
    InstantaneousZ.initialize();
  #endif
  ExploitMatrix.allocate(fyear,lyear,1,nages,"ExploitMatrix");
  #ifndef NO_AD_INITIALIZE
    ExploitMatrix.initialize();
  #endif
  PredAbundance.allocate(fyear,lyear+1,1,nages,"PredAbundance");
  #ifndef NO_AD_INITIALIZE
    PredAbundance.initialize();
  #endif
  PosfunPenalty.allocate("PosfunPenalty");
  #ifndef NO_AD_INITIALIZE
  PosfunPenalty.initialize();
  #endif
  FPosfunPenalty.allocate("FPosfunPenalty");
  #ifndef NO_AD_INITIALIZE
  FPosfunPenalty.initialize();
  #endif
  AgeVector.allocate(1,nages,"AgeVector");
  #ifndef NO_AD_INITIALIZE
    AgeVector.initialize();
  #endif
  PredKill.allocate(fyear,lyear,1,nages,"PredKill");
  #ifndef NO_AD_INITIALIZE
    PredKill.initialize();
  #endif
  KillSelectivity.allocate(1,nages,"KillSelectivity");
  #ifndef NO_AD_INITIALIZE
    KillSelectivity.initialize();
  #endif
  KillResiduals.allocate(fyear,lyear,1,nages,"KillResiduals");
  #ifndef NO_AD_INITIALIZE
    KillResiduals.initialize();
  #endif
  TotalKillResiduals.allocate(fyear,lyear,"TotalKillResiduals");
  #ifndef NO_AD_INITIALIZE
    TotalKillResiduals.initialize();
  #endif
  TotalKillLogLike.allocate("TotalKillLogLike");
  #ifndef NO_AD_INITIALIZE
  TotalKillLogLike.initialize();
  #endif
  KillAgeCompLogLike.allocate("KillAgeCompLogLike");
  #ifndef NO_AD_INITIALIZE
  KillAgeCompLogLike.initialize();
  #endif
  EstGillNetSelectivityParameters.allocate(1,4,"EstGillNetSelectivityParameters");
  #ifndef NO_AD_INITIALIZE
    EstGillNetSelectivityParameters.initialize();
  #endif
  GillNetSelectivity.allocate(1,nages,"GillNetSelectivity");
  #ifndef NO_AD_INITIALIZE
    GillNetSelectivity.initialize();
  #endif
  PredGillNet.allocate(fyear,lyear+1,1,nages,"PredGillNet");
  #ifndef NO_AD_INITIALIZE
    PredGillNet.initialize();
  #endif
  TotalGillNetResiduals.allocate(fyear,lyear+1,"TotalGillNetResiduals");
  #ifndef NO_AD_INITIALIZE
    TotalGillNetResiduals.initialize();
  #endif
  GillNetResiduals.allocate(fyear,lyear+1,1,nages,"GillNetResiduals");
  #ifndef NO_AD_INITIALIZE
    GillNetResiduals.initialize();
  #endif
  TotalGillNetLogLike.allocate("TotalGillNetLogLike");
  #ifndef NO_AD_INITIALIZE
  TotalGillNetLogLike.initialize();
  #endif
  GillNetAgeCompLogLike.allocate("GillNetAgeCompLogLike");
  #ifndef NO_AD_INITIALIZE
  GillNetAgeCompLogLike.initialize();
  #endif
  EstOffshoreGillNetSelectivityParameters.allocate(1,4,"EstOffshoreGillNetSelectivityParameters");
  #ifndef NO_AD_INITIALIZE
    EstOffshoreGillNetSelectivityParameters.initialize();
  #endif
  OffshoreGillNetSelectivity.allocate(1,nages,"OffshoreGillNetSelectivity");
  #ifndef NO_AD_INITIALIZE
    OffshoreGillNetSelectivity.initialize();
  #endif
  PredOffshoreGillNet.allocate(1999,lyear+1,1,nages,"PredOffshoreGillNet");
  #ifndef NO_AD_INITIALIZE
    PredOffshoreGillNet.initialize();
  #endif
  TotalOffshoreGillNetResiduals.allocate(1999,lyear+1,"TotalOffshoreGillNetResiduals");
  #ifndef NO_AD_INITIALIZE
    TotalOffshoreGillNetResiduals.initialize();
  #endif
  OffshoreGillNetResiduals.allocate(1999,lyear+1,1,nages,"OffshoreGillNetResiduals");
  #ifndef NO_AD_INITIALIZE
    OffshoreGillNetResiduals.initialize();
  #endif
  TotalOffshoreGillNetLogLike.allocate("TotalOffshoreGillNetLogLike");
  #ifndef NO_AD_INITIALIZE
  TotalOffshoreGillNetLogLike.initialize();
  #endif
  OffshoreGillNetAgeCompLogLike.allocate("OffshoreGillNetAgeCompLogLike");
  #ifndef NO_AD_INITIALIZE
  OffshoreGillNetAgeCompLogLike.initialize();
  #endif
  PredTrawlYOY.allocate(fyear,lyear+1,"PredTrawlYOY");
  #ifndef NO_AD_INITIALIZE
    PredTrawlYOY.initialize();
  #endif
  TrawlYOYWeights.allocate(fyear,lyear+1,"TrawlYOYWeights");
  #ifndef NO_AD_INITIALIZE
    TrawlYOYWeights.initialize();
  #endif
  TrawlYOYResiduals.allocate(fyear,lyear+1,"TrawlYOYResiduals");
  #ifndef NO_AD_INITIALIZE
    TrawlYOYResiduals.initialize();
  #endif
  TrawlYOYLogLike.allocate("TrawlYOYLogLike");
  #ifndef NO_AD_INITIALIZE
  TrawlYOYLogLike.initialize();
  #endif
  PredTrawlOnePlus.allocate(1989,lyear+1,"PredTrawlOnePlus");
  #ifndef NO_AD_INITIALIZE
    PredTrawlOnePlus.initialize();
  #endif
  TrawlOnePlusWeights.allocate(1989,lyear+1,"TrawlOnePlusWeights");
  #ifndef NO_AD_INITIALIZE
    TrawlOnePlusWeights.initialize();
  #endif
  TrawlOnePlusResiduals.allocate(1989,lyear+1,"TrawlOnePlusResiduals");
  #ifndef NO_AD_INITIALIZE
    TrawlOnePlusResiduals.initialize();
  #endif
  TrawlOnePlusLogLike.allocate("TrawlOnePlusLogLike");
  #ifndef NO_AD_INITIALIZE
  TrawlOnePlusLogLike.initialize();
  #endif
  SemiObsF.allocate(fyear,lyear,"SemiObsF");
  #ifndef NO_AD_INITIALIZE
    SemiObsF.initialize();
  #endif
  FResiduals.allocate(fyear,lyear,"FResiduals");
  #ifndef NO_AD_INITIALIZE
    FResiduals.initialize();
  #endif
  FLogLike.allocate("FLogLike");
  #ifndef NO_AD_INITIALIZE
  FLogLike.initialize();
  #endif
  LikelihoodLambda.allocate(1,9,"LikelihoodLambda");
  #ifndef NO_AD_INITIALIZE
    LikelihoodLambda.initialize();
  #endif
  NumberOfParameters.allocate("NumberOfParameters");
  #ifndef NO_AD_INITIALIZE
  NumberOfParameters.initialize();
  #endif
  AIC.allocate("AIC");
  #ifndef NO_AD_INITIALIZE
  AIC.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  getAnalysisData();
  getMissingDataIndicators();
  CalcObsAgeProportions();
  AgeVector.fill_seqadd(1,1);  // Vector with values 1 through nages to use in selectivity function calculations for gill net
}

void model_parameters::userfunction(void)
{
  f =0.0;
  NumberOfParameters=initial_params::nvarcalc();
  PosfunPenalty=0.0;
  FPosfunPenalty=0.0;
  TransformParameters();
  getKillSelectivities();
  getGillNetSelectivities();
  getOffshoreGillNetSelectivities();
  calcFZExp();
  calcPredAbundance();
  calcPlusAbundance();  //INVOKED TWICE INTENTIONALLY.  SEE IF ONE WILL WORK THE SAME.
  calcPredKill();
  calcPredGillNet();
  calcPredOffshoreGillNet();
  calcPredTrawl();
  CalcPredAgeProportions();  //WILL NEED TO THINK ABOUT ORDER OF FUNCTIONS AND OPERATIONS
  SemiObsF=Effortq*ObsEffort;
  TotalKillResiduals=log(ObsKillByYear)-log(PredKillByYear);
  TotalKillLogLike=-0.5*norm2(TotalKillResiduals/KillSigma);
  KillAgeCompLogLike=sum(KillSampSize*rowsum(elem_prod(ObsKillPropAge,log(PredKillPropAge))));
  KillResiduals=log(ObsKill)-log(PredKill);
  TotalGillNetResiduals=log(ObsGillNetByYear)-log(PredGillNetByYear);
  TotalGillNetLogLike=-0.5*norm2(TotalGillNetResiduals/GillNetSigma); 
  GillNetAgeCompLogLike=sum(GillNetSampSize*rowsum(elem_prod(ObsGillNetPropAge,log(PredGillNetPropAge))));
  GillNetResiduals=log(ObsGillNet)-log(PredGillNet);
  TotalOffshoreGillNetResiduals=log(ObsOffshoreGillNetByYear)-log(PredOffshoreGillNetByYear);
  TotalOffshoreGillNetLogLike=-0.5*norm2(TotalOffshoreGillNetResiduals/OffshoreGillNetSigma); 
  OffshoreGillNetAgeCompLogLike=sum(OffshoreGillNetSampSize*rowsum(elem_prod(ObsOffshoreGillNetPropAge,log(PredOffshoreGillNetPropAge))));
  OffshoreGillNetResiduals=log(ObsOffshoreGillNet)-log(PredOffshoreGillNet);
  TrawlYOYResiduals=log(ObsTrawlYOY)-log(PredTrawlYOY);
  TrawlYOYLogLike=-0.5*norm2(elem_prod(TrawlYOYResiduals,TrawlYOYWeights)/TrawlYOYSigma);
  TrawlOnePlusResiduals=log(ObsTrawlOnePlus)-log(PredTrawlOnePlus);
  TrawlOnePlusLogLike=-0.5*norm2(elem_prod(TrawlOnePlusResiduals,TrawlOnePlusWeights)/TrawlOnePlusSigma); 
  FResiduals=log(SemiObsF)-log(FByYear);
  FLogLike=-0.5*norm2(FResiduals/FSigma); //Divisor used as standard deviation
  calcPlusAbundance();
  LikelihoodLambda(1)=1;  //Lambda for total kill
  LikelihoodLambda(2)=1;  //Lambda for kill age compositions
  LikelihoodLambda(3)=1;  //Lambda for total gill net catch
  LikelihoodLambda(4)=1;  //Lambda for gill net age compositions
  LikelihoodLambda(5)=1;  //Lambda for total offshore gill net catch
  LikelihoodLambda(6)=1;  //Lambda for offshore gill net age compositions
  LikelihoodLambda(7)=1;  //Lambda for trawl YOY catch
  LikelihoodLambda(8)=1;  //Lambda for trawl age 1+ catch
  LikelihoodLambda(9)=1;  //Lambda for F - effort relationship
  f=-1.0*(TotalKillLogLike*LikelihoodLambda(1)
          +KillAgeCompLogLike*LikelihoodLambda(2)
          +TotalGillNetLogLike*LikelihoodLambda(3)
          +GillNetAgeCompLogLike*LikelihoodLambda(4)
          +TotalOffshoreGillNetLogLike*LikelihoodLambda(5)
          +OffshoreGillNetAgeCompLogLike*LikelihoodLambda(6)
          +TrawlYOYLogLike*LikelihoodLambda(7)
          +TrawlOnePlusLogLike*LikelihoodLambda(8)
          +FLogLike*LikelihoodLambda(9));  //Multiply entire quantity by -1 to find minimum of negative log likelihood, equivalent to max log likelihood
  AIC=(2*f)+(2*NumberOfParameters);
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::final_calcs()
{
  OutputFiles();
}

void model_parameters::getAnalysisData(void)
{
  ObsKill=ReadInObsKill.sub(fyear,lyear);
  ObsGillNet=ReadInObsGillNet.sub(fyear,lyear+1);
  ObsOffshoreGillNet=ReadInObsOffshoreGillNet.sub(1999,lyear+1);
  ObsTrawlYOY=ReadInObsTrawlYOY(fyear,lyear+1);
  ObsTrawlOnePlus=ReadInObsTrawlOnePlus(1989,lyear+1);
  ObsEffort=ReadInObsEffort(fyear,lyear);
}

void model_parameters::getMissingDataIndicators(void)
{
  int i;
  for (i=fyear;i<=lyear+1;i++)  //Code assumes data missing for entire years at a time, not just for one age
  {
    if (ObsTrawlYOY(i)!=99.99) {TrawlYOYWeights(i)=1;}
  }
  for (i=1989;i<=lyear+1;i++)  
  {
    if (ObsTrawlOnePlus(i)!=99.99) {TrawlOnePlusWeights(i)=1;}
  }
}

void model_parameters::TransformParameters(void)
{
  NAgeOne=mfexp(LogNAgeOne);
  Effortq=exp(LogEffortq);
  GillNetq=exp(LogGillNetq);
  OffshoreGillNetq=exp(LogOffshoreGillNetq);
  TrawlYOYq=exp(LogTrawlYOYq);
  TrawlOnePlusq=exp(LogTrawlOnePlusq);
}

void model_parameters::CalcObsAgeProportions(void)
{
  ObsKillByYear=rowsum(ObsKill);
  ObsGillNetByYear=rowsum(ObsGillNet);
  ObsOffshoreGillNetByYear=rowsum(ObsOffshoreGillNet);
  for (int i=fyear;i<=lyear+1;i++)
  {
    for (int j=1;j<=nages;j++)
    {
      if (i<=lyear) {ObsKillPropAge(i,j)=ObsKill(i,j)/ObsKillByYear(i);}
      ObsGillNetPropAge(i,j)=ObsGillNet(i,j)/ObsGillNetByYear(i);
      if (i>=1999) {ObsOffshoreGillNetPropAge(i,j)=ObsOffshoreGillNet(i,j)/ObsOffshoreGillNetByYear(i);}
    }
  }
}

void model_parameters::CalcPredAgeProportions(void)
{
  PredKillByYear=rowsum(PredKill);
  PredGillNetByYear=rowsum(PredGillNet);
  PredOffshoreGillNetByYear=rowsum(PredOffshoreGillNet);
  for (int i=fyear;i<=lyear+1;i++)
  {
    for (int j=1;j<=nages;j++)
    {
      if (i<=lyear) {PredKillPropAge(i,j)=PredKill(i,j)/PredKillByYear(i);}
      PredGillNetPropAge(i,j)=PredGillNet(i,j)/PredGillNetByYear(i);
      if (i>=1999) {PredOffshoreGillNetPropAge(i,j)=PredOffshoreGillNet(i,j)/PredOffshoreGillNetByYear(i);}
    }
  }
}

void model_parameters::getKillSelectivities(void)
{
  int KillSelectivityMethod;
  KillSelectivityMethod=3;  //Choose one of the options below
  if (KillSelectivityMethod==1) //  Option 1 for estimating a separate selectivity for all but one age
  {
    KillSelectivity(1)=EstKillSelectivity14(1);
    KillSelectivity(2)=EstKillSelectivity14(2);
    KillSelectivity(3)=EstKillSelectivity14(3);
    KillSelectivity(4)=EstKillSelectivity14(4);
    KillSelectivity(5)=1.0;
    KillSelectivity(6)=EstKillSelectivity6;
  }
  if (KillSelectivityMethod==3) //  Option 3 for three-parameter gamma function to estimate selectivities
    {
    MaxAgeMinusLocation=AgeVector-EstKillGammaSelectivityLocation;
      for (int i=1;i<=nages;i++)
      {
        if (MaxAgeMinusLocation(i)<=0.0) MaxAgeMinusLocation(i)=0.0000000001;
      }
      KillSelectivity=(1/exp(gammln(EstKillGammaSelectivityShape)))*(pow(EstKillGammaSelectivityScale,EstKillGammaSelectivityShape))*elem_prod(pow(MaxAgeMinusLocation,EstKillGammaSelectivityShape-1.0),
                               mfexp(-EstKillGammaSelectivityScale*(AgeVector-EstKillGammaSelectivityLocation)));
      KillSelectivity=KillSelectivity/KillSelectivity(5);
    }
}

void model_parameters::getGillNetSelectivities(void)
{
  int GillNetSelectivityMethod;
  GillNetSelectivityMethod=3;  //Choose one of the options below
  if (GillNetSelectivityMethod==1) //  Option 1 for estimating a separate selectivity for all but one age
  {
    GillNetSelectivity(1)=EstGillNetSelectivity13(1);
    GillNetSelectivity(2)=EstGillNetSelectivity13(2);
    GillNetSelectivity(3)=EstGillNetSelectivity13(3);
    GillNetSelectivity(4)=1.0;
    GillNetSelectivity(5)=EstGillNetSelectivity56(5);
    GillNetSelectivity(6)=EstGillNetSelectivity56(6);
  }
  if (GillNetSelectivityMethod==2) //  Option 2 for four-parameter double logistic function to estimate selectivities
  {
    EstGillNetSelectivityParameters(1)=EstGillNetSelectivityParameter1;
    EstGillNetSelectivityParameters(2)=EstGillNetSelectivityParameter2;
    EstGillNetSelectivityParameters(3)=EstGillNetSelectivityParameter3;
    EstGillNetSelectivityParameters(4)=EstGillNetSelectivityParameter4;
    GillNetSelectivity=elem_div((1-(1/(1+mfexp(EstGillNetSelectivityParameters(4)*(AgeVector-EstGillNetSelectivityParameters(3)))))),
                                (1+mfexp(-EstGillNetSelectivityParameters(2)*(AgeVector-EstGillNetSelectivityParameters(1)))));
    GillNetSelectivity=GillNetSelectivity/GillNetSelectivity(4);
  }
  if (GillNetSelectivityMethod==3) //  Option 3 for three-parameter gamma function to estimate selectivities
  {
      MaxAgeMinusLocation=AgeVector-EstGillNetGammaSelectivityLocation;
      for (int i=1;i<=nages;i++)
      {
        if (MaxAgeMinusLocation(i)<=0.0) MaxAgeMinusLocation(i)=0.0000000001;
      }
      GillNetSelectivity=(1/exp(gammln(EstGillNetGammaSelectivityShape)))*(pow(EstGillNetGammaSelectivityScale,EstGillNetGammaSelectivityShape))*elem_prod(pow(MaxAgeMinusLocation,EstGillNetGammaSelectivityShape-1.0),
                               mfexp(-EstGillNetGammaSelectivityScale*(AgeVector-EstGillNetGammaSelectivityLocation)));
      GillNetSelectivity=GillNetSelectivity/GillNetSelectivity(4);
  }
}

void model_parameters::getOffshoreGillNetSelectivities(void)
{
  int OffshoreGillNetSelectivityMethod;
  OffshoreGillNetSelectivityMethod=3;  //Choose one of the options below
  if (OffshoreGillNetSelectivityMethod==1)  // Option 1 for estimating a separate selectivity for all but one age
  {
    OffshoreGillNetSelectivity(1)=EstOffshoreGillNetSelectivity13(1);
    OffshoreGillNetSelectivity(2)=EstOffshoreGillNetSelectivity13(2);
    OffshoreGillNetSelectivity(3)=EstOffshoreGillNetSelectivity13(3);
    OffshoreGillNetSelectivity(4)=1.0;
    OffshoreGillNetSelectivity(5)=EstOffshoreGillNetSelectivity56(5);
    OffshoreGillNetSelectivity(6)=EstOffshoreGillNetSelectivity56(6);
  }
  if (OffshoreGillNetSelectivityMethod==2) // Option 2 for four-parameter double logistic function to estimate offshore gill net selectivities
  {
    EstOffshoreGillNetSelectivityParameters(1)=EstOffshoreGillNetSelectivityParameter1;
    EstOffshoreGillNetSelectivityParameters(2)=EstOffshoreGillNetSelectivityParameter2;
    EstOffshoreGillNetSelectivityParameters(3)=EstOffshoreGillNetSelectivityParameter3;
    EstOffshoreGillNetSelectivityParameters(4)=EstOffshoreGillNetSelectivityParameter4;
    OffshoreGillNetSelectivity=elem_div((1-(1/(1+mfexp(EstOffshoreGillNetSelectivityParameters(4)*(AgeVector-EstOffshoreGillNetSelectivityParameters(3)))))),
                                (1+mfexp(-EstOffshoreGillNetSelectivityParameters(2)*(AgeVector-EstOffshoreGillNetSelectivityParameters(1)))));
    OffshoreGillNetSelectivity=OffshoreGillNetSelectivity/OffshoreGillNetSelectivity(4);
  }
  if (OffshoreGillNetSelectivityMethod==3) //  Option 3 for three-parameter gamma function to estimate selectivities
  {
      MaxAgeMinusLocation=AgeVector-EstOffshoreGillNetGammaSelectivityLocation;
      for (int i=1;i<=nages;i++)
      {
        if (MaxAgeMinusLocation(i)<=0.0) MaxAgeMinusLocation(i)=0.0000000001;
      }
      OffshoreGillNetSelectivity=(1/exp(gammln(EstOffshoreGillNetGammaSelectivityShape)))*(pow(EstOffshoreGillNetGammaSelectivityScale,EstOffshoreGillNetGammaSelectivityShape))*elem_prod(pow(MaxAgeMinusLocation,EstOffshoreGillNetGammaSelectivityShape-1.0),
                               mfexp(-EstOffshoreGillNetGammaSelectivityScale*(AgeVector-EstOffshoreGillNetGammaSelectivityLocation)));
      OffshoreGillNetSelectivity=OffshoreGillNetSelectivity/OffshoreGillNetSelectivity(4);
  }
}

void model_parameters::calcPredAbundance(void)
{
  int i;
  int j;
  for (j=1;j<=nages;j++) // Row One 
  {
    PredAbundance(fyear,j)=NAgeOne(fyear+1-j);
  }
  for (i=fyear;i<=lyear+1;i++) // Column One
  {
    PredAbundance(i,1)=NAgeOne(i);
  }
  for (i=fyear+1;i<=lyear+1;i++) // Rest of table through one before last age
  {
   for (j=2;j<=nages-1;j++)
    {
      PredAbundance(i,j)=PredAbundance(i-1,j-1)*mfexp(-InstantaneousZ(i-1,j-1));
      PredAbundance(i,j)=posfun(PredAbundance(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }
  for (i=fyear+1;i<=lyear+1;i++) // Last column for last age
  {
    PredAbundance(i,nages)=(PredAbundance(i-1,nages-1)*mfexp(-InstantaneousZ(i-1,nages-1)))+
                   (PredAbundance(i-1,nages)*mfexp(-InstantaneousZ(i-1,nages)));
    PredAbundance(i,nages)=posfun(PredAbundance(i,nages),0.000001,PosfunPenalty);  //Constrain to be positive
  }
}

void model_parameters::calcFZExp(void)
{
  for (int i=fyear;i<=lyear;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      InstantaneousF(i,j)=FByYear(i)*KillSelectivity(j);  //DELETE LINE BELOW WHEN ALL IS WORKING RIGHT
      //InstantaneousF(i,j)=-log(posfun(((PredAbundance(i,j)-PredKill(i,j))*mfexp(-M(j)))/PredAbundance(i,j),0.000001,FPosfunPenalty))-M(j);
      InstantaneousZ(i,j)=InstantaneousF(i,j)+M(j);
      ExploitMatrix(i,j)=(InstantaneousF(i,j)/InstantaneousZ(i,j))*(1-exp(-InstantaneousZ(i,j)));
    }
  }
}

void model_parameters::calcPredKill(void)
{
  for (int i=fyear;i<=lyear;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      PredKill(i,j)=ExploitMatrix(i,j)*PredAbundance(i,j); //CAN PROBABLY DO WITHOUT LOOPS IF BELOW CONSTRAINT NOT NEEDED
      PredKill(i,j)=posfun(PredKill(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }
}

void model_parameters::calcPredGillNet(void)
{
  for (int i=fyear;i<=lyear+1;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      PredGillNet(i,j)=GillNetq*GillNetSelectivity(j)*PredAbundance(i,j);
      PredGillNet(i,j)=posfun(PredGillNet(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }
}

void model_parameters::calcPredOffshoreGillNet(void)
{
  for (int i=1999;i<=lyear+1;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      PredOffshoreGillNet(i,j)=OffshoreGillNetq*OffshoreGillNetSelectivity(j)*PredAbundance(i,j);
      PredOffshoreGillNet(i,j)=posfun(PredOffshoreGillNet(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }
}

void model_parameters::calcPredTrawl(void)
{
  for (int i=fyear;i<=lyear+1;i++) 
  {
    PredTrawlYOY(i)=TrawlYOYq*PredAbundance(i,1);
    PredTrawlYOY(i)=posfun(PredTrawlYOY(i),0.000001,PosfunPenalty);  //Constrain to be positive
    if (i>=1989){PredTrawlOnePlus(i)=TrawlOnePlusq*OnePlusAbundanceByYear(i);
    PredTrawlOnePlus(i)=posfun(PredTrawlOnePlus(i),0.000001,PosfunPenalty);}  //Constrain to be positive
  }
}

void model_parameters::calcPlusAbundance(void)
{
  int j;
  OnePlusAbundanceByYear=rowsum(PredAbundance);
  for (int i=fyear;i<=lyear+1;i++) 
  {
    ThreePlusAbundanceByYear(i)=0;
    for (j=3;j<=nages;j++)
    {
      ThreePlusAbundanceByYear(i)=PredAbundance(i,j)+ThreePlusAbundanceByYear(i);
    }
  }
}

void model_parameters::OutputFiles(void)
{
  ofstream outDataYears("c:\\DataAndWork\\Temp\\YepOutDataYears.txt",ios::trunc);
  outDataYears<<fyear<<endl;
  outDataYears<<lyear<<endl; //Do not include forecast year
  outDataYears<<1<<endl;
  outDataYears<<nages<<endl;
  outDataYears.close();
  ofstream outData1("c:\\DataAndWork\\Temp\\YepOutN.txt",ios::trunc);
  outData1<<PredAbundance<<endl;
  outData1.close(); 
  ofstream outData2("c:\\DataAndWork\\Temp\\YepOutTotKill.txt",ios::trunc);
  outData2<<ObsKillByYear<<endl;
  outData2<<PredKillByYear<<endl;
  outData2.close(); 
  ofstream outData3("c:\\DataAndWork\\Temp\\YepOutKillResid.txt",ios::trunc); 
  outData3<<KillResiduals<<endl; 
  outData3.close();
  ofstream outData4("c:\\DataAndWork\\Temp\\YepOutKillAgeProp.txt",ios::trunc); 
  outData4<<ObsKillPropAge<<endl; 
  outData4<<PredKillPropAge<<endl; 
  outData4.close();
  ofstream outData5("c:\\DataAndWork\\Temp\\YepOutKillSelect.txt",ios::trunc);
  outData5<<KillSelectivity<<endl;
  outData5.close();
  ofstream outData6("c:\\DataAndWork\\Temp\\YepOutTotGN.txt",ios::trunc);
  outData6<<ObsGillNetByYear<<endl;
  outData6<<PredGillNetByYear<<endl;
  outData6.close(); 
  ofstream outData7("c:\\DataAndWork\\Temp\\YepOutGNResid.txt",ios::trunc);
  outData7<<GillNetResiduals<<endl;
  outData7.close();
  ofstream outData8("c:\\DataAndWork\\Temp\\YepOutGNAgeProp.txt",ios::trunc); 
  outData8<<ObsGillNetPropAge<<endl; 
  outData8<<PredGillNetPropAge<<endl; 
  outData8.close();
  ofstream outData9("c:\\DataAndWork\\Temp\\YepOutGNSelect.txt",ios::trunc);
  outData9<<GillNetSelectivity<<endl;
  outData9.close();
  ofstream outData10("c:\\DataAndWork\\Temp\\YepOutTotOffGN.txt",ios::trunc);
  outData10<<ObsOffshoreGillNetByYear<<endl;
  outData10<<PredOffshoreGillNetByYear<<endl;
  outData10.close(); 
  ofstream outData11("c:\\DataAndWork\\Temp\\YepOutOffGNResid.txt",ios::trunc);
  outData11<<OffshoreGillNetResiduals<<endl;
  outData11.close();
  ofstream outData12("c:\\DataAndWork\\Temp\\YepOutOffGNAgeProp.txt",ios::trunc); 
  outData12<<ObsOffshoreGillNetPropAge<<endl; 
  outData12<<PredOffshoreGillNetPropAge<<endl; 
  outData12.close();
  ofstream outData13("c:\\DataAndWork\\Temp\\YepOutOffGNSelect.txt",ios::trunc);
  outData13<<OffshoreGillNetSelectivity<<endl;
  outData13.close();
  ofstream outData14("c:\\DataAndWork\\Temp\\YepOutTrawlYOY.txt",ios::trunc);
  outData14<<ObsTrawlYOY<<endl;
  outData14<<PredTrawlYOY<<endl;
  outData14.close();
  ofstream outData15("c:\\DataAndWork\\Temp\\YepOutTrawlYOYResid.txt",ios::trunc);
  outData15<<TrawlYOYResiduals<<endl;
  outData15.close();
  ofstream outData16("c:\\DataAndWork\\Temp\\YepOutTrawl1P.txt",ios::trunc);
  outData16<<ObsTrawlOnePlus<<endl;
  outData16<<PredTrawlOnePlus<<endl;
  outData16.close();
  ofstream outData17("c:\\DataAndWork\\Temp\\YepOutTrawl1PResid.txt",ios::trunc);
  outData17<<TrawlOnePlusResiduals<<endl;
  outData17.close();
  ofstream outData18("c:\\DataAndWork\\Temp\\YepOutFByYear.txt",ios::trunc);
  outData18<<SemiObsF<<endl;
  outData18<<FByYear<<endl;
  outData18.close(); 
  ofstream outData19("c:\\DataAndWork\\Temp\\YepOut3Plus.txt",ios::trunc);
  outData19<<ThreePlusAbundanceByYear<<endl;
  outData19.close(); 
  ofstream outData20("c:\\DataAndWork\\Temp\\YepOutInstF.txt",ios::trunc);
  outData20<<InstantaneousF<<endl;
  outData20.close();
  ofstream outData21("c:\\DataAndWork\\Temp\\YepOutEffort.txt",ios::trunc);
  outData21<<ObsEffort<<endl; 
  outData21<<FByYear/Effortq<<endl;
  outData21.close();
  ofstream outData22("c:\\DataAndWork\\Temp\\YepOutObsKill.txt",ios::trunc);
  outData22<<ObsKill<<endl; 
  outData22.close();
  ofstream outData23("c:\\DataAndWork\\Temp\\YepOutPredKill.txt",ios::trunc);
  outData23<<PredKill<<endl; 
  outData23.close();
  ofstream outData24("c:\\DataAndWork\\Temp\\YepOutObsGillNet.txt",ios::trunc);
  outData24<<ObsGillNet<<endl; 
  outData24.close();
  ofstream outData25("c:\\DataAndWork\\Temp\\YepOutPredGillNet.txt",ios::trunc);
  outData25<<PredGillNet<<endl; 
  outData25.close();
  ofstream outData26("c:\\DataAndWork\\Temp\\YepOutObsOffshoreGillNet.txt",ios::trunc);
  outData26<<ObsOffshoreGillNet<<endl; 
  outData26.close();
  ofstream outData27("c:\\DataAndWork\\Temp\\YepOutPredOffshoreGillNet.txt",ios::trunc);
  outData27<<PredOffshoreGillNet<<endl; 
  outData27.close();
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Predicted Abundance" << endl;
  report << PredAbundance << endl;
  report << "Three Plus Abundance By Year" << endl;
  report << ThreePlusAbundanceByYear << endl;
  report << "Instantaneous F Matrix" << endl;
  report << InstantaneousF << endl;
  report << "Observed Total Kill By Year" << endl;
  report << ObsKillByYear << endl;
  report << "Predicted Total Kill By Year" << endl;
  report << PredKillByYear << endl;
  report << "Total Kill By Year Log Residuals" << endl;
  report << TotalKillResiduals << endl;
  report << "Observed Kill By Year and Age" << endl;
  report << ObsKill << endl;
  report << "Predicted Kill By Year and Age" << endl;
  report << PredKill << endl;
  report << "Kill Log Residuals by Year and Age" << endl;
  report << KillResiduals << endl;
  report << "Observed Kill Proportion by Age" << endl;
  report << ObsKillPropAge << endl;
  report << "Predicted Kill Proportion by Age" << endl;
  report << PredKillPropAge << endl;
  report << "Kill Selectivities" << endl;
  report << KillSelectivity << endl;
  report << "Kill Gamma Selectivity Parameters (Shape, Scale, Location)" << endl;
  report << EstKillGammaSelectivityShape << " " << EstKillGammaSelectivityScale << " " << EstKillGammaSelectivityLocation<< endl;
  report << "Observed Total Gill Net Catch By Year" << endl;
  report << ObsGillNetByYear << endl;
  report << "Predicted Total Gill Net Catch By Year" << endl;
  report << PredGillNetByYear << endl;
  report << "Total Gill Net Catch By Year Log Residuals" << endl;
  report << TotalGillNetResiduals << endl;
  report << "Observed Gill Net by Year and Age" << endl;
  report << ObsGillNet << endl;
  report << "Predicted Gill Net by Year and Age" << endl;
  report << PredGillNet << endl;
  report << "Gill Net Log Residuals by Year and Age" << endl;
  report << GillNetResiduals << endl;
  report << "Observed Gill Net Proportion by Age" << endl;
  report << ObsGillNetPropAge << endl;
  report << "Predicted Gill Net Proportion by Age" << endl;
  report << PredGillNetPropAge << endl;
  report << "Gill Net Catchability Parameter q" << endl;
  report << GillNetq << endl;
  report << "Gill Net Selectivities" << endl;
  report << GillNetSelectivity << endl;
  report << "Gill Net Selectivity Parameters for Double Logistic Function" << endl;
  report << EstGillNetSelectivityParameters << endl;
  report << "Gill Net Gamma Selectivity Parameters (Shape, Scale, Location)" << endl;
  report << EstGillNetGammaSelectivityShape << " " << EstGillNetGammaSelectivityScale << " " << EstGillNetGammaSelectivityLocation<< endl;
  report << "Observed Total Offshore Gill Net Catch By Year" << endl;
  report << ObsOffshoreGillNetByYear << endl;
  report << "Predicted Total Offshore Gill Net Catch By Year" << endl;
  report << PredOffshoreGillNetByYear << endl;
  report << "Total Offshore Gill Net Catch By Year Log Residuals" << endl;
  report << TotalOffshoreGillNetResiduals << endl;
  report << "Observed Offshore Gill Net by Year and Age" << endl;
  report << ObsOffshoreGillNet << endl;
  report << "Predicted Offshore Gill Net by Year and Age" << endl;
  report << PredOffshoreGillNet << endl;
  report << "Offshore Gill Net Log Residuals by Year and Age" << endl;
  report << OffshoreGillNetResiduals << endl;
  report << "Observed Offshore Gill Net Proportion by Age" << endl;
  report << ObsOffshoreGillNetPropAge << endl;
  report << "Predicted Offshore Gill Net Proportion by Age" << endl;
  report << PredOffshoreGillNetPropAge << endl;
  report << "Offshore Gill Net Catchability Parameter q" << endl;
  report << OffshoreGillNetq << endl;
  report << "Offshore Gill Net Selectivities" << endl;
  report << OffshoreGillNetSelectivity << endl;
  report << "Offshore Gill Net Selectivity Parameters for Double Logistic Function" << endl;
  report << EstOffshoreGillNetSelectivityParameters << endl;
  report << "Offshore Gill Net Gamma Selectivity Parameters (Shape, Scale, Location)" << endl;
  report << EstOffshoreGillNetGammaSelectivityShape << " " << EstOffshoreGillNetGammaSelectivityScale << " " << EstOffshoreGillNetGammaSelectivityLocation<< endl;
  report << "TrawlYOYLogLike" << endl;
  report << TrawlYOYLogLike << endl;
  report << "TrawlOnePlusLogLike" << endl;
  report << TrawlOnePlusLogLike << endl;
  report << "Observed Trawl for YOY" << endl;
  report << ObsTrawlYOY << endl;
  report << "Predicted Trawl for YOY" << endl;
  report << PredTrawlYOY << endl;
  report << "Trawl YOY Residuals" << endl;
  report << TrawlYOYResiduals << endl;
  report << "Trawl YOY Catchability Parameter q" << endl;
  report << TrawlYOYq << endl;
  report << "Observed Trawl for Age 1+" << endl;
  report << ObsTrawlOnePlus << endl;
  report << "Predicted Trawl for Age 1+" << endl;
  report << PredTrawlOnePlus << endl;
  report << "Trawl Age 1+ Residuals" << endl;
  report << TrawlOnePlusResiduals << endl;
  report << "Trawl Age 1+ Catchability Parameter q" << endl;
  report << TrawlOnePlusq << endl;
  report << "Observed Effort" << endl;
  report << ObsEffort << endl;
  report << "Semi-Observed F (Observed Effort * Catchability q)" << endl;
  report << SemiObsF << endl;
  report << "Estimated Fishing Mortalities by Year" << endl;
  report << FByYear << endl;
  report << "F Log Residuals" << endl;
  report << FResiduals << endl;
  report << "PosfunPenalty (used to constrain abundances to be positive)" << endl;
  report << PosfunPenalty << endl;
  report << "FPosfunPenalty (used to prevent logs of negative numbers in calculating F" << endl;
  report << FPosfunPenalty << endl;
  report << "Total Kill Log Likelihood" << endl;
  report << TotalKillLogLike << endl;
  report << "Kill Age Composition Log Likelihood" << endl;
  report << KillAgeCompLogLike << endl;
  report << "Total Gill Net Catch Log Likelihood" << endl;
  report << TotalGillNetLogLike << endl;
  report << "Gill Net Age Composition Log Likelihood" << endl;
  report << GillNetAgeCompLogLike << endl;
  report << "Total Offshore Gill Net Catch Log Likelihood" << endl;
  report << TotalOffshoreGillNetLogLike << endl;
  report << "Offshore Gill Net Age Composition Log Likelihood" << endl;
  report << OffshoreGillNetAgeCompLogLike << endl;
  report << "Trawl YOY Catch Log Likelihood" << endl;
  report << TrawlYOYLogLike << endl;
  report << "Trawl Age 1+ Catch Log Likelihood" << endl;
  report << TrawlOnePlusLogLike << endl;
  report << "F Log Likelihood" << endl;
  report << FLogLike << endl;
  report << "Total Log Likelihood" << endl;
  report << -f << endl;
  report << "Number of Parameters" << endl;
  report << NumberOfParameters << endl;
  report << "AIC (the smaller the better)" << endl;
  report << AIC << endl;
  report << "Retrospective Analysis Quantities:  Age 3+, Age 1, Total Kill, F, Gill Net, Offshore Gill Net, Trawl YOY, Trawl 1+" << endl;
  report << ThreePlusAbundanceByYear << endl;
  report << NAgeOne(fyear,lyear+1) << endl;
  report << PredKillByYear << endl;
  report << FByYear << endl;
  report << PredGillNetByYear << endl;
  report << PredOffshoreGillNetByYear << endl;
  report << PredTrawlYOY << endl;
  report << PredTrawlOnePlus << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
