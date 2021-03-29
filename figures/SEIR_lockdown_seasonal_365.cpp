#include <Rcpp.h>
using namespace Rcpp;

// This contains the model functions for the SEIR model, with interaction
// in the I and R compartments

NumericVector transmission_calc_ER_SEIR_season (NumericVector y, NumericMatrix Contact_Structure, int num_grps, float bR)
{
  //initialise output vector
  NumericVector transmission_R(num_grps*22);
  //for the age group being infected
  for(int l=0;(num_grps)>l;++l){
    //set up as 0 initially
    transmission_R[l]=0;
    //for each infecting age group
    for (int k=0;(num_grps)>k;++k){
      //wotk out relative infection
      float temp = bR * ( y[2+k*22] * Contact_Structure(l, k) + y[6+k*22] * Contact_Structure(l, k) +
                            y[10+k*22] * Contact_Structure(l, k) + y[14+k*22] * Contact_Structure(l, k));
      //add to previous infections of age group
      transmission_R[l]= transmission_R[l]+temp;
    };
  }; return(transmission_R);
}

NumericVector transmission_calc_EI_SEIR_season (NumericVector y, NumericMatrix Contact_Structure, int num_grps, float bI)
{
  //initialise output vector
  NumericVector transmission_I(num_grps*22);
  //for the age group being infected
  for(int a=0;(num_grps)>a;++a){
    //set up as 0 initially
    //  Rprintf("%i a", a);
    transmission_I[a]=0;
    //Rprintf("\n here 2");
    //for each infecting age group
    for (int b=0;(num_grps)>b;++b){
      //wotk out relative infection
      float temp = bI * ( y[8+b*22] * Contact_Structure(a, b) + y[9+b*22] * Contact_Structure(a, b) +
                            y[10+b*22] * Contact_Structure(a, b) + y[11+b*22] * Contact_Structure(a, b));
      //add to previous infections of age group
      transmission_I[a]= transmission_I[a]+temp;
      // Rprintf("\n transmission I is %f and temp is %f",transmission_I[a],temp);
    };
  }; return(transmission_I);
}

//seasonal forcing function
float seasonal_other_SEIR_season (double year_time, float amplitude, float baseline,
                                  float phi, float gamma)
{
  float seasonal_beta = (amplitude* std::cos(2*M_PI*(year_time-phi)/365.25) +
                           baseline);
  return(seasonal_beta);
}

//[[Rcpp::export]]
List SEIR_2virus_cons_season (double t, NumericVector y, List parms)
{
  
  // parameters
  int num_grps = parms["num_grps"];
  float bR = parms["beta_covid_0"];
  // float bR_0 = parms["beta_covid_0"];
  // float bR_1 = parms["beta_covid_1"];
  float bI = parms["beta_other"];
  float sigma_R = parms["sig_R"];
  float sigma_I = parms["sig_I"];
  float sigma_solid_R = parms["sig_solid_I"];// note, these were specified the wrong way around so intentially wswapping
  float sigma_solid_I = parms["sig_solid_R"];
  //  Rcout << sigma_I << "\n";
  float gammaR =  parms["gamma_covid"];
  float gammaI = parms["gamma_other"];
  // float rho =  parms["rho"];
  
  // NumericVector trickleI_t = parms["other_intros"];
  NumericVector trickleR_t = parms["covid_intros"];
  float waningR = parms["waning_covid"];
  float waningI = parms["waning_other"];
  //  float N = parms["N"];
  float incubation_R = parms["incubation_covid"];
  float incubation_I = parms["incubation_other"];
  float reporting_delay_R1 = parms["inf_to_symp_covid"];
  float reporting_delay_I1 = parms["inf_to_symp_other"];
  float reporting_delay_R2 = parms["reporting_delay_covid"];
  float reporting_delay_I2 = parms["reporting_delay_other"];
  
  NumericVector RSV_Sus = parms["red_susc_covid"];
  NumericMatrix contacts_hh = parms["contacts_hh"];
  NumericMatrix contacts_other = parms["contacts_other"];
  NumericMatrix contacts_school = parms["contacts_school"];
  NumericVector ydot(num_grps*22);
  float covid_time = parms["covid_time"];
  float amplitude = parms["amplitude"];
  float amplitude_covid = parms["amplitude_covid"];
  float phi = parms["phi"];
  float mobility_start = parms["mobility_start"];
  // float child_extra_reduc = parms["child_extra_reduc"];
  
  // ageing rates
  NumericVector ageing_rates = parms["age_rates"];
  float births = parms["births"];
  float deaths = parms["deaths"];
  
  // float bR;
  float R0change1 = parms["covid_change_time"];
  
  NumericVector trickleR (num_grps);
  NumericMatrix Contact_Structure (num_grps) ;
  NumericMatrix c_all = parms["contacts_all"] ;
  // NumericMatrix Contact_Structure = c_all;
  NumericVector contacts_reduc = parms["contacts_reduc"];
  float social_distancing = parms["social_distancing"];
  float start_dying_off = parms["start_dying_off"];
  // NumericVector trickleI;
  
  if(t<mobility_start){
    //hh, school and other contacts
    for (int te = 0; te < 16; te++) {
      for(int to = 0; to < 16; to++){
        Contact_Structure(te,to) =  c_all(te,to);
      }}
  } else if(t>=mobility_start){
    for (int te = 0; te < 16; te++) {
      for(int to = 0; to < 16; to++){
        //hh contacts + mobility other + schools
        Contact_Structure(te,to) =  c_all(te,to) -
          //          proprortion of all contacts in ohter category * %change in contacts * actual numbers
        (contacts_other(te,to)*(1-contacts_reduc[t-mobility_start])*(c_all(te,to)));
      }}
    // schools closed on top
    if(t>= (R0change1 + covid_time)){
      for (int te = 0; te < 16; te++) {
        for(int to = 0; to < 16; to++){
          //  the alerady caluclated reduction from mobiliuty - the proportion at school*total contacts
          Contact_Structure(te,to) = Contact_Structure(te,to) - (contacts_school(te,to)*c_all(te,to));
        }}
      
      Contact_Structure = Contact_Structure * social_distancing;
    }}
  
  // vector of 0s
  NumericVector trickleR2 (num_grps);
  
  if ((t > covid_time) & (t < (covid_time + R0change1))){
    trickleR = trickleR_t;
    //      Rcout << trickleR[1] << "\n";
  } else{
    trickleR = trickleR2;
  }
  
  bI = seasonal_other_SEIR_season(t,amplitude, bI, phi, gammaI);
  bR = seasonal_other_SEIR_season(t,amplitude_covid, bR, phi, gammaR);
  
  // calculate the population size.
  float N = 0;
  float counters = 0;
  /// sum all the counting states that are not included in the total population number
  for( int q=0; (num_grps)>q; ++q){
    counters = counters + y[16+q*22];
    counters = counters + y[17+q*22];
    counters = counters + y[18+q*22];
    counters = counters + y[19+q*22];
    counters = counters + y[20+q*22];
    counters = counters + y[21+q*22];
  }
  
  // calculate total number of people at dying age
  float dying_age = 0;
  for( int r=start_dying_off; (num_grps)>r; ++r){
    
    for(int sub_r=0; 16>sub_r; ++sub_r){
      dying_age = dying_age + y[sub_r+r*22];
    }}
  
  //Rcout << dying_age << "\n";
  N = sum(y) - counters ;
  //  Rcout << N << "\n";
  NumericVector transmission_R = transmission_calc_ER_SEIR_season (y, Contact_Structure,num_grps,bR);
  NumericVector transmission_I = transmission_calc_EI_SEIR_season (y, Contact_Structure,num_grps,bI);
  
  for( int i=0; (num_grps)>i; ++i){
    
    // for the youngest age group
    
    if(i==0){
      
      //SS0
      ydot[0+i*22]= ( - (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      - (transmission_I[i]*y[0+i*22])
                      + waningI*y[12+i*22]
                      + waningR*y[3+i*22]
                      - trickleR[i]
                      - ageing_rates[i]* y[0+i*22]
                      + births
      );
      
      //ES
      ydot[1+i*22]= ( + (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      + waningI*y[13+i*22]
                      - (incubation_R*y[1+i*22])
                      - (transmission_I[i]*y[1+i*22])
                      - ageing_rates[i]* y[1+i*22]
                      + trickleR[i]
      );
      
      //IS
      ydot[2+i*22] = ( + ((1-sigma_solid_R)*incubation_R*y[1+i*22])
                       + waningI*y[14+i*22]
                       - (gammaR*y[2+i*22])
                       - (sigma_I*transmission_I[i]*y[2+i*22])
                       - ageing_rates[i]* y[2+i*22]
      );
      
      //RS
      ydot[3+i*22] = ( + (gammaR*y[2+i*22])
                       + waningI*y[15+i*22]
                       - waningR*y[3+i*22]
                       - (sigma_I*transmission_I[i]*y[3+i*22])
                       - ageing_rates[i]* y[3+i*22]
      );
      
      
      //SE
      ydot[4+i*22] = ( + waningR* y[7+i*22]
                       + transmission_I[i]*y[0+i*22]
                       - (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       - (incubation_I *y[4+i*22])
                       - ageing_rates[i]* y[4+i*22]
                       //       + trickleI[i]
      );
      
      
      //EE
      ydot[5+i*22] = ( + (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       + transmission_I[i]*y[1+i*22]
                       - (incubation_R*y[5+i*22])
                       - (incubation_I*y[5+i*22])
                       - ageing_rates[i]* y[5+i*22]
      );
      
      //IE
      ydot[6+i*22] = ( + incubation_R*y[5+i*22]
                       + sigma_I*transmission_I[i]*y[2+i*22]
                       - (gammaR*y[6+i*22])
                       - (incubation_I*y[6+i*22])
                       - ageing_rates[i]* y[6+i*22]
      );
      
      //RE
      ydot[7+i*22] = ( + gammaR*y[6+i*22]
                       + sigma_I*transmission_I[i]*y[3+i*22]
                       - waningR*y[7+i*22]
                       - incubation_I*y[7+i*22]
                       - ageing_rates[i]* y[7+i*22]
      );
      
      
      //SI
      ydot[8+i*22] = ( + waningR*y[11+i*22]
                       + (1-sigma_solid_I)*incubation_I*y[4+i*22]
                       - (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                       - (gammaI*y[8+i*22])
                       - ageing_rates[i]* y[8+i*22]
      );
      
      
      //EI
      ydot[9+i*22] = ( + (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                       + incubation_I*y[5+i*22]
                       - incubation_R*y[9+i*22]
                       - gammaI*y[9+i*22]
                       - ageing_rates[i]* y[9+i*22]
      );
      
      //II
      ydot[10+i*22] = (+ incubation_R*y[9+i*22]
                       + incubation_I*y[6+i*22]
                       - (gammaR*y[10+i*22])
                       - (gammaI*y[10+i*22])
                       - ageing_rates[i]* y[10+i*22]
      );
      
      //RI
      ydot[11+i*22] = (+ gammaR*y[10+i*22]
                       + incubation_I*y[7+i*22]
                       - waningR*y[11+i*22]
                       - gammaI*y[11+i*22]
                       - ageing_rates[i]* y[11+i*22]
                       + sigma_solid_I*incubation_I*(y[4+i*22])
      );
      
      //SR
      ydot[12+i*22] = ( + waningR*y[15+i*22]
                        + gammaI*y[8+i*22]
                        - (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                        - waningI*y[12+i*22]
                        - ageing_rates[i]* y[12+i*22]
      );
      
      
      //ER
      ydot[13+i*22] = ( + (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                        + gammaR*y[9+i*22]
                        - incubation_R*y[13+i*22]
                        - waningI*y[13+i*22]
                        - ageing_rates[i]* y[13+i*22]
      );
      
      //IR
      ydot[14+i*22] = ( + incubation_R*y[13+i*22]
                        + gammaI*y[10+i*22]
                        - gammaR*y[14+i*22]
                        - waningI*y[14+i*22]
                        - ageing_rates[i]* y[14+i*22]
                        + (sigma_solid_R*incubation_R*(y[1+i*22]))
      );
      
      //RR
      ydot[15+i*22] = ( + gammaI*y[11+i*22]
                        + gammaR*y[14+i*22]
                        - waningR*y[15+i*22]
                        - waningI*y[15+i*22]
                        - ageing_rates[i]* y[15+i*22]
      );
      
      
      //R cases
      ydot[16+i*22] = (+ incubation_R*y[1+i*22]
                       + incubation_R*y[5+i*22]
                       + incubation_R*y[9+i*22]
                       + incubation_R*y[13+i*22]
                       - y[16+i*22]*reporting_delay_R1
      );
      
      
      //I cases
      ydot[17+i*22] =( + incubation_I*y[4+i*22]
                       + incubation_I*y[5+i*22]
                       + incubation_I*y[6+i*22]
                       + incubation_I*y[7+i*22]
                       - y[17+i*22]*reporting_delay_I1
      );
      
      ydot[18+i*22] = y[16+i*22]*reporting_delay_R1
      -y[18+i*22]*reporting_delay_R2;
      ydot[19+i*22] = y[17+i*22]*reporting_delay_I1
      -y[19+i*22]*reporting_delay_I2;
      
      ydot[20+i*22] = y[18+i*22]*reporting_delay_R2;
      ydot[21+i*22] = y[19+i*22]*reporting_delay_I2;
      
      
    } else if (i==(num_grps - 1)){
      // for the last age group - only dying out
      //SS0
      ydot[0+i*22]= ( - (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      - (transmission_I[i]*y[0+i*22])
                      + waningI*y[12+i*22]
                      + waningR*y[3+i*22]
                      - trickleR[i]
                      + ageing_rates[i - 1]* y[0+(i - 1)*22]
                      - deaths*y[0+i*22]/dying_age
      );
      
      //ES
      ydot[1+i*22]= ( + (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      + waningI*y[13+i*22]
                      - (incubation_R*y[1+i*22])
                      - (transmission_I[i]*y[1+i*22])
                      + ageing_rates[i - 1]* y[1+(i - 1)*22]
                      - deaths*y[1+i*22]/dying_age
                      + trickleR[i]
      );
      
      //IS
      ydot[2+i*22] = ( + ((1-sigma_solid_R)*incubation_R*y[1+i*22])
                       + waningI*y[14+i*22]
                       - (gammaR*y[2+i*22])
                       - (sigma_I*transmission_I[i]*y[2+i*22])
                       + ageing_rates[i - 1]* y[2+(i - 1)*22]
                       - deaths*y[2+i*22]/dying_age
      );
      
      //RS
      ydot[3+i*22] = ( + (gammaR*y[2+i*22])
                       + waningI*y[15+i*22]
                       - waningR*y[3+i*22]
                       - (sigma_I*transmission_I[i]*y[3+i*22])
                       + ageing_rates[i - 1]* y[3+(i - 1)*22]
                       - deaths*y[3+i*22]/dying_age
      );
      
      
      //SE
      ydot[4+i*22] = ( + waningR* y[7+i*22]
                       + transmission_I[i]*y[0+i*22]
                       - (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       - (incubation_I *y[4+i*22])
                       + ageing_rates[i - 1]* y[4+(i - 1)*22]
                       - deaths*y[4+i*22]/dying_age
                       //       + trickleI[i]
      );
      
      
      //EE
      ydot[5+i*22] = ( + (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       + transmission_I[i]*y[1+i*22]
                       - (incubation_R*y[5+i*22])
                       - (incubation_I*y[5+i*22])
                       + ageing_rates[i - 1]* y[5+(i - 1)*22]
                       - deaths*y[5+i*22]/dying_age
      );
      
      //IE
      ydot[6+i*22] = ( + incubation_R*y[5+i*22]
                       + sigma_I*transmission_I[i]*y[2+i*22]
                       - (gammaR*y[6+i*22])
                       - (incubation_I*y[6+i*22])
                       + ageing_rates[i - 1]* y[6+(i - 1)*22]
                       - deaths*y[6+i*22]/dying_age
      );
      
      //RE
      ydot[7+i*22] = ( + gammaR*y[6+i*22]
                       + sigma_I*transmission_I[i]*y[3+i*22]
                       - waningR*y[7+i*22]
                       - incubation_I*y[7+i*22]
                       + ageing_rates[i - 1]* y[7+(i - 1)*22]
                       - deaths*y[7+i*22]/dying_age
      );
      
      
      //SI
      ydot[8+i*22] = ( + waningR*y[11+i*22]
                       + (1-sigma_solid_I)*incubation_I*y[4+i*22]
                       - (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                       - (gammaI*y[8+i*22])
                       + ageing_rates[i - 1]* y[8+(i - 1)*22]
                       - deaths*y[8+i*22]/dying_age
      );
      
      
      //EI
      ydot[9+i*22] = (  + (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                        + incubation_I*y[5+i*22]
                        - incubation_R*y[9+i*22]
                        - gammaI*y[9+i*22]
                        + ageing_rates[i - 1]* y[9+(i - 1)*22]
                        - deaths*y[9+i*22]/dying_age
      );
      
      //II
      ydot[10+i*22] = (+ incubation_R*y[9+i*22]
                       + incubation_I*y[6+i*22]
                       - (gammaR*y[10+i*22])
                       - (gammaI*y[10+i*22])
                       + ageing_rates[i - 1]* y[10+(i - 1)*22]
                       - deaths*y[10+i*22]/dying_age
      );
      
      //RI
      ydot[11+i*22] = (+ gammaR*y[10+i*22]
                       + incubation_I*y[7+i*22]
                       - waningR*y[11+i*22]
                       - gammaI*y[11+i*22]
                       + ageing_rates[i - 1]* y[11+(i - 1)*22]
                       - deaths*y[11+i*22]/dying_age
                       + sigma_solid_I*incubation_I*(y[4+i*22])
                       
      );
      
      //SR
      ydot[12+i*22] = ( + waningR*y[15+i*22]
                        + gammaI*y[8+i*22]
                        - (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                        - waningI*y[12+i*22]
                        + ageing_rates[i - 1]* y[12+(i - 1)*22]
                        - deaths*y[12+i*22]/dying_age
      );
      
      
      //ER
      ydot[13+i*22] = (  + (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                         + gammaR*y[9+i*22]
                         - incubation_R*y[13+i*22]
                         - waningI*y[13+i*22]
                         + ageing_rates[i - 1]* y[13+(i - 1)*22]
                         - deaths*y[13+i*22]/dying_age
      );
      
      //IR
      ydot[14+i*22] = ( + incubation_R*y[13+i*22]
                        + gammaI*y[10+i*22]
                        - gammaR*y[14+i*22]
                        - waningI*y[14+i*22]
                        + ageing_rates[i - 1]* y[14+(i - 1)*22]
                        - deaths*y[14+i*22]/dying_age
                        + (sigma_solid_R*incubation_R*(y[1+i*22]))
                        
      );
      
      //RR
      ydot[15+i*22] = ( + gammaI*y[11+i*22]
                        + gammaR*y[14+i*22]
                        - waningR*y[15+i*22]
                        - waningI*y[15+i*22]
                        + ageing_rates[i - 1]* y[15+(i - 1)*22]
                        - deaths*y[15+i*22]/dying_age
      );
      
      
      //R cases
      ydot[16+i*22] = (+ incubation_R*y[1+i*22]
                       + incubation_R*y[5+i*22]
                       + incubation_R*y[9+i*22]
                       + incubation_R*y[13+i*22]
                       - y[16+i*22]*reporting_delay_R1
      );
      
      
      //I cases
      ydot[17+i*22] =( + incubation_I*y[4+i*22]
                       + incubation_I*y[5+i*22]
                       + incubation_I*y[6+i*22]
                       + incubation_I*y[7+i*22]
                       - y[17+i*22]*reporting_delay_I1
      );
      
      ydot[18+i*22] = y[16+i*22]*reporting_delay_R1
      -y[18+i*22]*reporting_delay_R2;
      ydot[19+i*22] = y[17+i*22]*reporting_delay_I1
      -y[19+i*22]*reporting_delay_I2;
      
      ydot[20+i*22] = y[18+i*22]*reporting_delay_R2;
      ydot[21+i*22] = y[19+i*22]*reporting_delay_I2;
      
      // all age groups above 60 - have aging in, aging out and dying
    } else if (i>=start_dying_off){
      
      //SS0
      ydot[0+i*22]= ( - (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      - (transmission_I[i]*y[0+i*22])
                      + waningI*y[12+i*22]
                      + waningR*y[3+i*22]
                      - trickleR[i]
                      + ageing_rates[i - 1]* y[0+(i - 1)*22]
                      - ageing_rates[i]*y[0 + i*22]
                      - deaths*y[0+i*22]/dying_age
      );
      
      //ES
      ydot[1+i*22]= ( + (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      + waningI*y[13+i*22]
                      - (incubation_R*y[1+i*22])
                      - (transmission_I[i]*y[1+i*22])
                      + ageing_rates[i - 1]* y[1+(i - 1)*22]
                      + trickleR[i]
                      - ageing_rates[i]*y[1 + i*22]
                      - deaths*y[1+i*22]/dying_age
      );
      
      //IS
      ydot[2+i*22] = ( + ((1-sigma_solid_R)*incubation_R*y[1+i*22])
                       + waningI*y[14+i*22]
                       - (gammaR*y[2+i*22])
                       - (sigma_I*transmission_I[i]*y[2+i*22])
                       + ageing_rates[i - 1]* y[2+(i - 1)*22]
                       - ageing_rates[i]*y[2 + i*22]
                       - deaths*y[2+i*22]/dying_age
      );
      
      //RS
      ydot[3+i*22] = ( + (gammaR*y[2+i*22])
                       + waningI*y[15+i*22]
                       - waningR*y[3+i*22]
                       - (sigma_I*transmission_I[i]*y[3+i*22])
                       + ageing_rates[i - 1]* y[3+(i - 1)*22]
                       - ageing_rates[i]*y[3 + i*22]
                       - deaths*y[3+i*22]/dying_age
      );
      
      
      //SE
      ydot[4+i*22] = ( + waningR* y[7+i*22]
                       + transmission_I[i]*y[0+i*22]
                       - (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       - (incubation_I *y[4+i*22])
                       + ageing_rates[i - 1]* y[4+(i - 1)*22]
                       - ageing_rates[i]*y[4 + i*22]
                       - deaths*y[4+i*22]/dying_age
                       //       + trickleI[i]
      );
      
      
      //EE
      ydot[5+i*22] = ( + (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       + transmission_I[i]*y[1+i*22]
                       - (incubation_R*y[5+i*22])
                       - (incubation_I*y[5+i*22])
                       + ageing_rates[i - 1]* y[5+(i - 1)*22]
                       - ageing_rates[i]*y[5 + i*22]
                       - deaths*y[5+i*22]/dying_age
      );
      
      //IE
      ydot[6+i*22] = ( + incubation_R*y[5+i*22]
                       + sigma_I*transmission_I[i]*y[2+i*22]
                       - (gammaR*y[6+i*22])
                       - (incubation_I*y[6+i*22])
                       + ageing_rates[i - 1]* y[6+(i - 1)*22]
                       - ageing_rates[i]*y[6 + i*22]
                       - deaths*y[6+i*22]/dying_age
      );
      
      //RE
      ydot[7+i*22] = ( + gammaR*y[6+i*22]
                       + sigma_I*transmission_I[i]*y[3+i*22]
                       - waningR*y[7+i*22]
                       - incubation_I*y[7+i*22]
                       + ageing_rates[i - 1]* y[7+(i - 1)*22]
                       - ageing_rates[i]*y[7 + i*22]
                       - deaths*y[7+i*22]/dying_age
      );
      
      
      //SI
      ydot[8+i*22] = ( + waningR*y[11+i*22]
                       + (1-sigma_solid_I)*incubation_I*y[4+i*22]
                       - (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                       - (gammaI*y[8+i*22])
                       + ageing_rates[i - 1]* y[8+(i - 1)*22]
                       - ageing_rates[i]*y[8 + i*22]
                       - deaths*y[8+i*22]/dying_age
      );
      
      
      //EI
      ydot[9+i*22] = ( + (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                       + incubation_I*y[5+i*22]
                       - incubation_R*y[9+i*22]
                       - gammaI*y[9+i*22]
                       + ageing_rates[i - 1]* y[9+(i - 1)*22]
                       - ageing_rates[i]*y[9 + i*22]
                       - deaths*y[9+i*22]/dying_age
      );
      
      //II
      ydot[10+i*22] = (+ incubation_R*y[9+i*22]
                       + incubation_I*y[6+i*22]
                       - (gammaR*y[10+i*22])
                       - (gammaI*y[10+i*22])
                       + ageing_rates[i - 1]* y[10+(i - 1)*22]
                       - ageing_rates[i]*y[10 + i*22]
                       - deaths*y[10+i*22]/dying_age
      );
      
      //RI
      ydot[11+i*22] = (+ gammaR*y[10+i*22]
                       + incubation_I*y[7+i*22]
                       - waningR*y[11+i*22]
                       - gammaI*y[11+i*22]
                       + ageing_rates[i - 1]* y[11+(i - 1)*22]
                       - ageing_rates[i]*y[11 + i*22]
                       - deaths*y[11+i*22]/dying_age
                       + sigma_solid_I*incubation_I*(y[4+i*22])
                       
      );
      
      //SR
      ydot[12+i*22] = ( + waningR*y[15+i*22]
                        + gammaI*y[8+i*22]
                        - (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                        - waningI*y[12+i*22]
                        + ageing_rates[i - 1]* y[12+(i - 1)*22]
                        - ageing_rates[i]*y[12 + i*22]
                        - deaths*y[12+i*22]/dying_age
      );
      
      
      //ER
      ydot[13+i*22] = (  + (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                         + gammaR*y[9+i*22]
                         - incubation_R*y[13+i*22]
                         - waningI*y[13+i*22]
                         + ageing_rates[i - 1]* y[13+(i - 1)*22]
                         - ageing_rates[i]*y[13 + i*22]
                         - deaths*y[13+i*22]/dying_age
      );
      
      //IR
      ydot[14+i*22] = ( + incubation_R*y[13+i*22]
                        + gammaI*y[10+i*22]
                        - gammaR*y[14+i*22]
                        - waningI*y[14+i*22]
                        + ageing_rates[i - 1]* y[14+(i - 1)*22]
                        - ageing_rates[i]*y[14 + i*22]
                        - deaths*y[14+i*22]/dying_age
                        + (sigma_solid_R*incubation_R*(y[1+i*22]))
                        
      );
      
      //RR
      ydot[15+i*22] = ( + gammaI*y[11+i*22]
                        + gammaR*y[14+i*22]
                        - waningR*y[15+i*22]
                        - waningI*y[15+i*22]
                        + ageing_rates[i - 1]* y[15+(i - 1)*22]
                        - ageing_rates[i]*y[15 + i*22]
                        - deaths*y[15+i*22]/dying_age
      );
      
      
      //R cases
      ydot[16+i*22] = (+ incubation_R*y[1+i*22]
                       + incubation_R*y[5+i*22]
                       + incubation_R*y[9+i*22]
                       + incubation_R*y[13+i*22]
                       - y[16+i*22]*reporting_delay_R1
      );
      
      
      //I cases
      ydot[17+i*22] =( + incubation_I*y[4+i*22]
                       + incubation_I*y[5+i*22]
                       + incubation_I*y[6+i*22]
                       + incubation_I*y[7+i*22]
                       - y[17+i*22]*reporting_delay_I1
      );
      
      ydot[18+i*22] = y[16+i*22]*reporting_delay_R1
      -y[18+i*22]*reporting_delay_R2;
      ydot[19+i*22] = y[17+i*22]*reporting_delay_I1
      -y[19+i*22]*reporting_delay_I2;
      
      ydot[20+i*22] = y[18+i*22]*reporting_delay_R2;
      ydot[21+i*22] = y[19+i*22]*reporting_delay_I2;
      
    }else{
      //all in between age groups
      //SS0
      ydot[0+i*22]= ( - (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      - (transmission_I[i]*y[0+i*22])
                      + waningI*y[12+i*22]
                      + waningR*y[3+i*22]
                      - trickleR[i]
                      + ageing_rates[i - 1]* y[0+(i - 1)*22]
                      - ageing_rates[i]*y[0 + i*22]
      );
      
      //ES
      ydot[1+i*22]= ( + (RSV_Sus[i]*transmission_R[i]*y[0+i*22])
                      + waningI*y[13+i*22]
                      - (incubation_R*y[1+i*22])
                      - (transmission_I[i]*y[1+i*22])
                      + ageing_rates[i - 1]* y[1+(i - 1)*22]
                      + trickleR[i]
                      - ageing_rates[i]*y[1 + i*22]
      );
      
      //IS
      ydot[2+i*22] = ( + ((1-sigma_solid_R)*incubation_R*y[1+i*22])
                       + waningI*y[14+i*22]
                       - (gammaR*y[2+i*22])
                       - (sigma_I*transmission_I[i]*y[2+i*22])
                       + ageing_rates[i - 1]* y[2+(i - 1)*22]
                       - ageing_rates[i]*y[2 + i*22]
      );
      
      //RS
      ydot[3+i*22] = ( + (gammaR*y[2+i*22])
                       + waningI*y[15+i*22]
                       - waningR*y[3+i*22]
                       - (sigma_I*transmission_I[i]*y[3+i*22])
                       + ageing_rates[i - 1]* y[3+(i - 1)*22]
                       - ageing_rates[i]*y[3 + i*22]
      );
      
      
      //SE
      ydot[4+i*22] = ( + waningR* y[7+i*22]
                       + transmission_I[i]*y[0+i*22]
                       - (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       - (incubation_I *y[4+i*22])
                       + ageing_rates[i - 1]* y[4+(i - 1)*22]
                       - ageing_rates[i]*y[4 + i*22]
                       //       + trickleI[i]
      );
      
      
      //EE
      ydot[5+i*22] = ( + (RSV_Sus[i]*transmission_R[i]*y[4+i*22])
                       + transmission_I[i]*y[1+i*22]
                       - (incubation_R*y[5+i*22])
                       - (incubation_I*y[5+i*22])
                       + ageing_rates[i - 1]* y[5+(i - 1)*22]
                       - ageing_rates[i]*y[5 + i*22]
      );
      
      //IE
      ydot[6+i*22] = (  + incubation_R*y[5+i*22]
                        + sigma_I*transmission_I[i]*y[2+i*22]
                        - (gammaR*y[6+i*22])
                        - (incubation_I*y[6+i*22])
                        + ageing_rates[i - 1]* y[6+(i - 1)*22]
                        - ageing_rates[i]*y[6 + i*22]
      );
      
      //RE
      ydot[7+i*22] = ( + gammaR*y[6+i*22]
                       + sigma_I*transmission_I[i]*y[3+i*22]
                       - waningR*y[7+i*22]
                       - incubation_I*y[7+i*22]
                       + ageing_rates[i - 1]* y[7+(i - 1)*22]
                       - ageing_rates[i]*y[7 + i*22]
      );
      
      
      //SI
      ydot[8+i*22] = ( + waningR*y[11+i*22]
                       + (1-sigma_solid_I)*incubation_I*y[4+i*22]
                       - (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                       - (gammaI*y[8+i*22])
                       + ageing_rates[i - 1]* y[8+(i - 1)*22]
                       - ageing_rates[i]*y[8 + i*22]
      );
      
      
      //EI
      ydot[9+i*22] = (  + (sigma_R*RSV_Sus[i]*transmission_R[i]*(y[8+i*22]))
                        + incubation_I*y[5+i*22]
                        - incubation_R*y[9+i*22]
                        - gammaI*y[9+i*22]
                        + ageing_rates[i - 1]* y[9+(i - 1)*22]
                        - ageing_rates[i]*y[9 + i*22]
      );
      
      //II
      ydot[10+i*22] = (+ incubation_R*y[9+i*22]
                       + incubation_I*y[6+i*22]
                       - (gammaR*y[10+i*22])
                       - (gammaI*y[10+i*22])
                       + ageing_rates[i - 1]* y[10+(i - 1)*22]
                       - ageing_rates[i]*y[10 + i*22]
      );
      
      //RI
      ydot[11+i*22] = (+ gammaR*y[10+i*22]
                       + incubation_I*y[7+i*22]
                       - waningR*y[11+i*22]
                       - gammaI*y[11+i*22]
                       + ageing_rates[i - 1]* y[11+(i - 1)*22]
                       - ageing_rates[i]*y[11 + i*22]
                       + sigma_solid_I*incubation_I*(y[4+i*22])
                       
      );
      
      //SR
      ydot[12+i*22] = ( + waningR*y[15+i*22]
                        + gammaI*y[8+i*22]
                        - (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                        - waningI*y[12+i*22]
                        + ageing_rates[i - 1]* y[12+(i - 1)*22]
                        - ageing_rates[i]*y[12 + i*22]
      );
      
      
      //ER
      ydot[13+i*22] = (  + (sigma_R*RSV_Sus[i]*transmission_R[i]*y[12+i*22])
                         + gammaR*y[9+i*22]
                         - incubation_R*y[13+i*22]
                         - waningI*y[13+i*22]
                         + ageing_rates[i - 1]* y[13+(i - 1)*22]
                         - ageing_rates[i]*y[13 + i*22]
      );
      
      //IR
      ydot[14+i*22] = ( + incubation_R*y[13+i*22]
                        + gammaI*y[10+i*22]
                        - gammaR*y[14+i*22]
                        - waningI*y[14+i*22]
                        + ageing_rates[i - 1]* y[14+(i - 1)*22]
                        - ageing_rates[i]*y[14 + i*22]
                        + (sigma_solid_R*incubation_R*(y[1+i*22]))
                        
      );
      
      //RR
      ydot[15+i*22] = ( + gammaI*y[11+i*22]
                        + gammaR*y[14+i*22]
                        - waningR*y[15+i*22]
                        - waningI*y[15+i*22]
                        + ageing_rates[i - 1]* y[15+(i - 1)*22]
                        - ageing_rates[i]*y[15 + i*22]
      );
      
      
      //R cases
      ydot[16+i*22] = (+ incubation_R*y[1+i*22]
                       + incubation_R*y[5+i*22]
                       + incubation_R*y[9+i*22]
                       + incubation_R*y[13+i*22]
                       - y[16+i*22]*reporting_delay_R1
      );
      
      
      //I cases
      ydot[17+i*22] =( + incubation_I*y[4+i*22]
                       + incubation_I*y[5+i*22]
                       + incubation_I*y[6+i*22]
                       + incubation_I*y[7+i*22]
                       - y[17+i*22]*reporting_delay_I1
      );
      
      ydot[18+i*22] = y[16+i*22]*reporting_delay_R1
      -y[18+i*22]*reporting_delay_R2;
      ydot[19+i*22] = y[17+i*22]*reporting_delay_I1
      -y[19+i*22]*reporting_delay_I2;
      
      ydot[20+i*22] = y[18+i*22]*reporting_delay_R2;
      ydot[21+i*22] = y[19+i*22]*reporting_delay_I2;
      
    }
    
  }
  
  return List::create(_["yout"] = ydot) ;
}

