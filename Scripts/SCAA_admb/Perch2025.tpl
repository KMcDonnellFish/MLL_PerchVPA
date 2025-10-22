// Mille Lacs Lake yellow perch statistical catch at age model, Rick Madsen, GLIFWC, August 2016

DATA_SECTION
  init_int fyear //Year indexes the year that the cohort is age 1
  init_int lyearread //Last year of harvest data as read in
  init_int lyear //Last year of harvest data for analysis. lyear+1 is year of forecast.
  init_int nages
  init_matrix ReadInObsKill(fyear,lyearread,1,nages) //Kill years include open water harvest for the year and for the preceding winter
  init_matrix ReadInObsGillNet(fyear,lyearread+1,1,nages) //Indices for survey years are incremented by 1 to refer to year of prediction.  Ages not incremented.
  init_matrix ReadInObsOffshoreGillNet(1999,lyearread+1,1,nages)
  init_vector ReadInObsTrawlYOY(fyear,lyearread+2)       //Trawl YOY years incremented by 2 to refer to year of age 1 harvest
  init_vector ReadInObsTrawlOnePlus(1989,lyearread+1)    //Trawl Age 1+ years incremented by 1 like gill net survey data to refer to year of prediction.
  init_vector M(1,nages)
  init_vector ReadInObsEffort(fyear,lyearread)
  init_number KillSigma
  init_number GillNetSigma
  init_number OffshoreGillNetSigma
  init_number TrawlYOYSigma
  init_number TrawlOnePlusSigma
  init_number FSigma
  init_number KillSampSize
  init_number GillNetSampSize
  init_number OffshoreGillNetSampSize

INITIALIZATION_SECTION
  LogNAgeOne 14.5
  FByYear 0.5
  EstKillGammaSelectivityShape 4.0
  EstKillGammaSelectivityScale 0.7
  EstGillNetGammaSelectivityShape 2.0
  EstGillNetGammaSelectivityScale 0.2
  EstOffshoreGillNetGammaSelectivityShape 2.0
  EstOffshoreGillNetGammaSelectivityScale 0.2

PARAMETER_SECTION  //Phases set to -1 if parameter not estimated
  init_bounded_vector LogNAgeOne(fyear-nages+1,lyear+1,7.6,20,1) //Abundance at age 1 (or first age) parameters

  init_bounded_vector FByYear(fyear,lyear,0.01,2,2) //Estimated fishing effort on age of reference selectivity
  init_bounded_number LogEffortq(-25,-2,1) //Effort-F catchability constant

  init_bounded_vector EstKillSelectivity14(1,4,0,5,-3) //Estimated kill selectivities for ages 1 to 4
  init_bounded_number EstKillSelectivity6(0,5,-3) //Estimated kill selectivities for age 6+
  init_bounded_number EstKillGammaSelectivityShape(1.1,10,3) //Shape parameter for gamma selectivity function
  init_bounded_number EstKillGammaSelectivityScale(0.001,5,3) //Scale parameter for gamma selectivity function
  init_bounded_number EstKillGammaSelectivityLocation(-1.0,1.0,-3) //Location parameter for gamma selectivity function

  init_bounded_number LogGillNetq(-25,-2,1) //Assessment gill net catchability
  init_bounded_vector EstGillNetSelectivity13(1,3,0,5,-3) //Estimated assessment gill net selectivities for ages 1 to 3
  init_bounded_vector EstGillNetSelectivity56(5,6,0,5,-3) //Estimated assessment gill net selectivities for ages 5 and 6+
  init_bounded_number EstGillNetSelectivityParameter1(-5,10,-3) // For double-logistic selectivity function, age of inflection for first curve
  init_bounded_number EstGillNetSelectivityParameter2(-2.5,2.5,-3) // For double-logistic selectivity function, slope for first curve
  init_bounded_number EstGillNetSelectivityParameter3(-5,10,-3) // For double-logistic selectivity function, age of inflection for second curve
  init_bounded_number EstGillNetSelectivityParameter4(-2.5,2.5,-3) // For double-logistic selectivity function, slope for second curve
  init_bounded_number EstGillNetGammaSelectivityShape(1.1,10,3) //Shape parameter for gamma selectivity function
  init_bounded_number EstGillNetGammaSelectivityScale(0.001,5,3) //Scale parameter for gamma selectivity function
  init_bounded_number EstGillNetGammaSelectivityLocation(-1.0,1.0,-3) //Location parameter for gamma selectivity function

  init_bounded_number LogOffshoreGillNetq(-25,-2,1) //Offshore assessment gill net catchability
  init_bounded_vector EstOffshoreGillNetSelectivity13(1,3,0,5,-3) //Estimated offshore assessment gill net selectivities for ages 1 to 3
  init_bounded_vector EstOffshoreGillNetSelectivity56(5,6,0,5,-3) //Estimated offshore assessment gill net selectivities for ages 5 and 6+
  init_bounded_number EstOffshoreGillNetSelectivityParameter1(-5,10,-3) // For double-logistic selectivity function, age of inflection for first curve
  init_bounded_number EstOffshoreGillNetSelectivityParameter2(-2.5,2.5,-3) // For double-logistic selectivity function, slope for first curve
  init_bounded_number EstOffshoreGillNetSelectivityParameter3(-5,10,-3) // For double-logistic selectivity function, age of inflection for second curve
  init_bounded_number EstOffshoreGillNetSelectivityParameter4(-2.5,2.5,-3) // For double-logistic selectivity function, slope for second curve
  init_bounded_number EstOffshoreGillNetGammaSelectivityShape(1.1,10,3) //Shape parameter for gamma selectivity function
  init_bounded_number EstOffshoreGillNetGammaSelectivityScale(0.001,5,3) //Scale parameter for gamma selectivity function
  init_bounded_number EstOffshoreGillNetGammaSelectivityLocation(-1.0,1.0,-3) //Location parameter for gamma selectivity function

  init_bounded_number LogTrawlYOYq(-25,-2,1) //Assessment trawl catchability for age 0 data
  init_bounded_number LogTrawlOnePlusq(-25,-2,1) //Assessment trawl catchability for age 1+ data

  vector NAgeOne(fyear-nages+1,lyear+1) //For transforming parameters estimated on log scale
  number Effortq
  number GillNetq
  number OffshoreGillNetq
  number TrawlYOYq
  number TrawlOnePlusq

  matrix ObsKill(fyear,lyear,1,nages)  //For getting subsets of data for retrospective analyses
  matrix ObsGillNet(fyear,lyear+1,1,nages)
  matrix ObsOffshoreGillNet(1999,lyear+1,1,nages)
  vector ObsTrawlYOY(fyear,lyear+1)
  vector ObsTrawlOnePlus(1989,lyear+1)
  vector ObsEffort(fyear,lyear)

  vector ObsKillByYear(fyear,lyear)  //For observed totals by year
  vector ObsGillNetByYear(fyear,lyear+1)
  vector ObsOffshoreGillNetByYear(1999,lyear+1)

  matrix ObsKillPropAge(fyear,lyear,1,nages)  //For observed proportions at age
  matrix ObsGillNetPropAge(fyear,lyear+1,1,nages)
  matrix ObsOffshoreGillNetPropAge(1999,lyear+1,1,nages)

  vector PredKillByYear(fyear,lyear)  //For predicted totals by year
  vector PredGillNetByYear(fyear,lyear+1)
  vector PredOffshoreGillNetByYear(1999,lyear+1)

  matrix PredKillPropAge(fyear,lyear,1,nages)  //For predicted proportions at age
  matrix PredGillNetPropAge(fyear,lyear+1,1,nages)
  matrix PredOffshoreGillNetPropAge(1999,lyear+1,1,nages)

  vector OnePlusAbundanceByYear(fyear,lyear+1)
  vector ThreePlusAbundanceByYear(fyear,lyear+1)
  vector MaxAgeMinusLocation(1,nages) // Needed for gamma calculations

  matrix InstantaneousF(fyear,lyear,1,nages) //F matrix
  matrix InstantaneousZ(fyear,lyear,1,nages) //Z matrix
  matrix ExploitMatrix(fyear,lyear,1,nages) //Exploitation rate matrix

  matrix PredAbundance(fyear,lyear+1,1,nages) //Predicted abundance matrix
  number PosfunPenalty //Cumulative "penalty" when posfun function invoked to constrain abundance to be positive
  number FPosfunPenalty //Cumulative "penalty" when posfun function used in calculation of instantaneous fishing mortalities

  vector AgeVector(1,nages) //Vector of ages to use in selectivity calculations
  matrix PredKill(fyear,lyear,1,nages) //Predicted kill matrix
  vector KillSelectivity(1,nages) //Kill selectivities for all ages
  matrix KillResiduals(fyear,lyear,1,nages) //Kill residuals for year and age
  vector TotalKillResiduals(fyear,lyear) //Total kill residuals by year
  number TotalKillLogLike //Kill log likelihood for total catch
  number KillAgeCompLogLike //Kill log likelihood for age compositions

  vector EstGillNetSelectivityParameters(1,4) //Vector to hold double-logistic function parameters
  vector GillNetSelectivity(1,nages) //Assessment gill net selectivities for all ages
  matrix PredGillNet(fyear,lyear+1,1,nages) //Predicted assessment gill net matrix
  vector TotalGillNetResiduals(fyear,lyear+1) //Total gill net residuals by year
  matrix GillNetResiduals(fyear,lyear+1,1,nages) //Assessment gill net residuals
  number TotalGillNetLogLike //Assessment gill net log likelihood for total catch
  number GillNetAgeCompLogLike //Assessment gill net log likelihood for age compositions

  vector EstOffshoreGillNetSelectivityParameters(1,4) //Vector to hold double-logistic function parameters
  vector OffshoreGillNetSelectivity(1,nages) //Assessment gill net selectivities for all ages
  matrix PredOffshoreGillNet(1999,lyear+1,1,nages) //Predicted assessment gill net matrix
  vector TotalOffshoreGillNetResiduals(1999,lyear+1) //Total gill net residuals by year
  matrix OffshoreGillNetResiduals(1999,lyear+1,1,nages) //Assessment gill net residuals
  number TotalOffshoreGillNetLogLike //Assessment gill net log likelihood for total catch
  number OffshoreGillNetAgeCompLogLike //Assessment gill net log likelihood for age compositions

  vector PredTrawlYOY(fyear,lyear+1) //Predicted trawl YOY (indexed here as age 1)
  vector TrawlYOYWeights(fyear,lyear+1) //Weights to account for missing data
  vector TrawlYOYResiduals(fyear,lyear+1) //Trawl YOY residuals
  number TrawlYOYLogLike //Trawl YOY log likelihood
  vector PredTrawlOnePlus(1989,lyear+1) //Predicted trawl for age 1+ (indexed here as age 2+)
  vector TrawlOnePlusWeights(1989,lyear+1) //Weights to account for missing data
  vector TrawlOnePlusResiduals(1989,lyear+1) //Trawl age 1+ residuals
  number TrawlOnePlusLogLike //Trawl age 1+ log likelihood

  vector SemiObsF(fyear,lyear) //Predicted effort
  vector FResiduals(fyear,lyear) //Effort residuals
  number FLogLike //Log likelihood for effort - F relationship

  vector LikelihoodLambda(1,9) //Vector to store lambdas for likelihood components
  number NumberOfParameters //Number of parameters in model
  number AIC //Akaike's information criterion
  objective_function_value f

PRELIMINARY_CALCS_SECTION

  getAnalysisData();
  getMissingDataIndicators();
  CalcObsAgeProportions();

  AgeVector.fill_seqadd(1,1);  // Vector with values 1 through nages to use in selectivity function calculations for gill net

PROCEDURE_SECTION

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

FINAL_SECTION

  OutputFiles();

FUNCTION getAnalysisData
  ObsKill=ReadInObsKill.sub(fyear,lyear);
  ObsGillNet=ReadInObsGillNet.sub(fyear,lyear+1);
  ObsOffshoreGillNet=ReadInObsOffshoreGillNet.sub(1999,lyear+1);
  ObsTrawlYOY=ReadInObsTrawlYOY(fyear,lyear+1);
  ObsTrawlOnePlus=ReadInObsTrawlOnePlus(1989,lyear+1);
  ObsEffort=ReadInObsEffort(fyear,lyear);

FUNCTION getMissingDataIndicators
  int i;
  for (i=fyear;i<=lyear+1;i++)  //Code assumes data missing for entire years at a time, not just for one age
  {
    if (ObsTrawlYOY(i)!=99.99) {TrawlYOYWeights(i)=1;}
  }
  for (i=1989;i<=lyear+1;i++)  
  {
    if (ObsTrawlOnePlus(i)!=99.99) {TrawlOnePlusWeights(i)=1;}
  }

FUNCTION TransformParameters  //For transforming parameters estimated on log scale back to regular
  NAgeOne=mfexp(LogNAgeOne);
  Effortq=exp(LogEffortq);
  GillNetq=exp(LogGillNetq);
  OffshoreGillNetq=exp(LogOffshoreGillNetq);
  TrawlYOYq=exp(LogTrawlYOYq);
  TrawlOnePlusq=exp(LogTrawlOnePlusq);

FUNCTION CalcObsAgeProportions
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

FUNCTION CalcPredAgeProportions
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

FUNCTION getKillSelectivities
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

FUNCTION getGillNetSelectivities
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

FUNCTION getOffshoreGillNetSelectivities
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

FUNCTION calcPredAbundance
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

FUNCTION calcFZExp
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

FUNCTION calcPredKill
  for (int i=fyear;i<=lyear;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      PredKill(i,j)=ExploitMatrix(i,j)*PredAbundance(i,j); //CAN PROBABLY DO WITHOUT LOOPS IF BELOW CONSTRAINT NOT NEEDED
      PredKill(i,j)=posfun(PredKill(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }

FUNCTION calcPredGillNet
  for (int i=fyear;i<=lyear+1;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      PredGillNet(i,j)=GillNetq*GillNetSelectivity(j)*PredAbundance(i,j);
      PredGillNet(i,j)=posfun(PredGillNet(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }

FUNCTION calcPredOffshoreGillNet
  for (int i=1999;i<=lyear+1;i++) 
  {
    for (int j=1;j<=nages;j++)
    {
      PredOffshoreGillNet(i,j)=OffshoreGillNetq*OffshoreGillNetSelectivity(j)*PredAbundance(i,j);
      PredOffshoreGillNet(i,j)=posfun(PredOffshoreGillNet(i,j),0.000001,PosfunPenalty);  //Constrain to be positive
    }
  }

FUNCTION calcPredTrawl
  for (int i=fyear;i<=lyear+1;i++) 
  {
    PredTrawlYOY(i)=TrawlYOYq*PredAbundance(i,1);
    PredTrawlYOY(i)=posfun(PredTrawlYOY(i),0.000001,PosfunPenalty);  //Constrain to be positive
    if (i>=1989){PredTrawlOnePlus(i)=TrawlOnePlusq*OnePlusAbundanceByYear(i);
    PredTrawlOnePlus(i)=posfun(PredTrawlOnePlus(i),0.000001,PosfunPenalty);}  //Constrain to be positive
  }

FUNCTION calcPlusAbundance  //Get age 1+ and 3+ summaries
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

FUNCTION OutputFiles
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

REPORT_SECTION
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
