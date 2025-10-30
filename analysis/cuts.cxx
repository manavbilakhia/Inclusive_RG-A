
#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TRandom3.h>
#include <fstream>
#include <stdio.h>


// IT IS IMPORTANT TO SET IT EVERY TIME DF APPLIES THE CUTS
// IT IS USED FOR SF CUT ONLY (SO FAR)

// bool isMC = false;

//Root file and plot path 
//string filePath = "files/";
//string plotPath = "files/";

//Mass and Energy parameters
const double m_p  = 0.93827208816;          //proton mass
const double b_E  = 10.6041;                //beam energy fall 2018
const double m_e  = 0.5109989461 * 0.001;   //electron mass
const double m_pi = 0.134977;               //pi0 mass

/////////////////// get 4-momentum, apply mom corr if it is data: ///////////////////////

inline auto dppC(float Px, float Py, float Pz, int sec, int ivec, int corEl, int corPip, int corPim, int corPro){
    
    // 'Px'/'Py'/'Pz'   ==> Corresponds to the Cartesian Components of the particle momentum being corrected
    // 'sec'            ==> Corresponds to the Forward Detector Sectors where the given particle is detected (6 total)
    // 'ivec'           ==> Corresponds to the particle being corrected (See below)    
        // (*) ivec = 0 --> Electron Corrections
        // (*) ivec = 1 --> Pi+ Corrections
        // (*) ivec = 2 --> Pi- Corrections
        // (*) ivec = 3 --> Proton Corrections
    // 'corEl'/'corPip'/'corPim'/'corPro' ==> Controls which version of the particle correction is used
        // Includes:
            // (*) Correction On/Off
            // (*) Pass Version
            // (*) Data Set (Fall 2018 or Spring 2019)
    // 'corEl'         ==> Controls the ELECTRON Corrections
        // corEl == 0  --> No Correction (Off)
        // corEl == 1  --> Fall  2018 - Pass 1
        // corEl == 2  --> Sping 2019 - Pass 2
        // corEl == 3  --> Fall  2018 - Pass 2
    // 'corPip'        ==> Controls the π+ PION Corrections
        // corPip == 0 --> No Correction
        // corPip == 1 --> Fall  2018 - Pass 1
        // corPip == 2 --> Sping 2019 - Pass 2
        // corPip == 3 --> Fall  2018 - Pass 2
    // 'corPim'        ==> Controls the π- PION Corrections
        // corPim == 0 --> No Correction
        // corPim == 1 --> Fall  2018 - Pass 1 (Created by Nick Trotta)
    // 'corPro'        ==> Controls the PROTON Corrections (Momentum)
        // corPro == 0 --> No Correction
        // corPro == 1 --> Fall  2018 - Pass 1

    // Momentum Magnitude
    double pp = sqrt(Px*Px + Py*Py + Pz*Pz);

    // Initializing the correction factor
    double dp = 0;

    // Defining Phi Angle
    double Phi = (180/3.1415926)*atan2(Py, Px);

    // Central Detector Corrections Not Included (Yet)

    // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
    if(((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)){
        Phi += 360;
    }

    // Getting Local Phi Angle
    double PhiLocal = Phi - (sec - 1)*60;

    // Applying Shift Functions to Phi Angles (local shifted phi = phi)
    double phi = PhiLocal;

    // For Electron Shift
    if(ivec == 0){
        phi = PhiLocal - 30/pp;
    }

    // For π+ Pion/Proton Shift
    if(ivec == 1 || ivec == 3){
        phi = PhiLocal + (32/(pp-0.05));
    }

    // For π- Pion Shift
    if(ivec == 2){
        phi = PhiLocal - (32/(pp-0.05));
    }


    //===============//===============//     No Corrections     //===============//===============//
    if(corEl == 0 && ivec == 0){ // No Electron Correction
        return dp/pp;
    }
    if(corPip == 0 && ivec == 1){ // No π+ Pion Correction
        return dp/pp;
    }
    if(corPim == 0 && ivec == 2){ // No π- Pion Correction
        return dp/pp;
    }
    if(corPro == 0 && ivec == 3){ // No Proton Correction
        return dp/pp;
    }
    //==============//==============//     No Corrections (End)     //==============//==============//

    //==============================//     Electron Corrections     //==============================//
    if(corEl != 0 && ivec == 0){
        if(corEl == 1){ // Fall 2018 - Pass 1 Corrections
            if(sec == 1){
                dp = ((-4.3303e-06)*phi*phi +  (1.1006e-04)*phi + (-5.7235e-04))*pp*pp +  ((3.2555e-05)*phi*phi +  (-0.0014559)*phi +   (0.0014878))*pp + ((-1.9577e-05)*phi*phi +   (0.0017996)*phi + (0.025963));
            }
            if(sec == 2){
                dp = ((-9.8045e-07)*phi*phi +  (6.7395e-05)*phi + (-4.6757e-05))*pp*pp + ((-1.4958e-05)*phi*phi +  (-0.0011191)*phi +  (-0.0025143))*pp +  ((1.2699e-04)*phi*phi +   (0.0033121)*phi + (0.020819));
            }
            if(sec == 3){
                dp = ((-5.9459e-07)*phi*phi + (-2.8289e-05)*phi + (-4.3541e-04))*pp*pp + ((-1.5025e-05)*phi*phi +  (5.7730e-04)*phi +  (-0.0077582))*pp +  ((7.3348e-05)*phi*phi +   (-0.001102)*phi + (0.057052));
            }
            if(sec == 4){
                dp = ((-2.2714e-06)*phi*phi + (-3.0360e-05)*phi + (-8.9322e-04))*pp*pp +  ((2.9737e-05)*phi*phi +  (5.1142e-04)*phi +   (0.0045641))*pp + ((-1.0582e-04)*phi*phi + (-5.6852e-04)*phi + (0.027506));
            }
            if(sec == 5){
                dp = ((-1.1490e-06)*phi*phi + (-6.2147e-06)*phi + (-4.7235e-04))*pp*pp +  ((3.7039e-06)*phi*phi + (-1.5943e-04)*phi + (-8.5238e-04))*pp +  ((4.4069e-05)*phi*phi +   (0.0014152)*phi + (0.031933));
            }
            if(sec == 6){
                dp =  ((1.1076e-06)*phi*phi +  (4.0156e-05)*phi + (-1.6341e-04))*pp*pp + ((-2.8613e-05)*phi*phi + (-5.1861e-04)*phi +  (-0.0056437))*pp +  ((1.2419e-04)*phi*phi +  (4.9084e-04)*phi + (0.049976));
            }
        }

        if(corEl == 2){ // Spring 2019 - Pass 2 Corrections
            if(sec == 1){
                dp = ((-1.4215599999999998e-06)*phi*phi + (4.91084e-06)*phi + (-0.00012995999999999998))*pp*pp + ((1.6952059999999994e-05)*phi*phi + (-0.00033224299999999997)*phi + (-0.0018080400000000003))*pp + ((-3.1853499999999996e-05)*phi*phi + (0.0016001439999999997)*phi + (0.03187985));
            }
            if(sec == 2){
                dp =             ((-5.4471e-06)*phi*phi + (-4.69579e-05)*phi + (0.000462807))*pp*pp + ((5.0258819999999995e-05)*phi*phi + (0.00023192399999999994)*phi + (-0.01118006))*pp + ((-8.5754e-05)*phi*phi + (0.00017097299999999994)*phi + (0.05324023));
            }
            if(sec == 3){
                dp = ((-3.4392460000000005e-06)*phi*phi + (9.860100000000002e-06)*phi + (-3.8414000000000015e-05))*pp*pp + ((1.7492300000000002e-05)*phi*phi + (-4.111499999999996e-05)*phi + (-0.0052975509999999984))*pp + ((1.0045499999999984e-05)*phi*phi + (-2.1412000000000004e-05)*phi + (0.03514576));
            }
            if(sec == 4){
                dp =  ((2.4865599999999998e-06)*phi*phi + (2.9090599999999996e-05)*phi + (0.00016154500000000003))*pp*pp + ((-2.7148730000000002e-05)*phi*phi + (-0.000136352)*phi + (-0.00543832))*pp + ((4.917660000000001e-05)*phi*phi + (-0.0001558459999999999)*phi + (0.04322285000000001));
            }
            if(sec == 5){
                dp =  ((-5.340280000000001e-06)*phi*phi + (1.355319e-05)*phi + (0.001362661))*pp*pp + ((5.858976999999999e-05)*phi*phi + (-0.00024119909999999995)*phi + (-0.02025752))*pp + ((-0.0001475504)*phi*phi + (0.0005707250000000001)*phi + (0.07970399));
            }
            if(sec == 6){
                dp = ((-3.0325500000000003e-06)*phi*phi + (-4.7810870999999994e-05)*phi + (0.001092504))*pp*pp + ((2.4123071999999996e-05)*phi*phi + (0.00047091400000000007)*phi + (-0.01504266))*pp + ((-9.523899999999999e-06)*phi*phi + (-0.0008819019999999998)*phi + (0.048088700000000005));
            }
        }
        
        if(corEl == 3){ // Fall 2018 - Pass 2 Corrections
            if(sec == 1){
                dp            =                ((-9.82416e-06)*phi*phi +            (-2.29956e-05)*phi +  (0.00029664199999999996))*pp*pp +           ((0.0001113414)*phi*phi +  (-2.041300000000001e-05)*phi +            (-0.00862226))*pp +            ((-0.000281738)*phi*phi +             (0.00058712)*phi +              (0.0652737));
                if(pp < 7){dp = dp +            ((-3.4001e-06)*phi*phi +             (-2.2885e-05)*phi +              (9.9705e-04))*pp*pp +             ((2.1840e-05)*phi*phi +              (2.4238e-04)*phi +             (-0.0091904))*pp +             ((-2.9180e-05)*phi*phi +            (-6.4496e-04)*phi +               (0.022505));}
                else{      dp = dp +            ((-6.3656e-05)*phi*phi +              (1.7266e-04)*phi +              (-0.0017909))*pp*pp +                ((0.00104)*phi*phi +              (-0.0028401)*phi +                (0.02981))*pp +              ((-0.0041995)*phi*phi +               (0.011537)*phi +                (-0.1196));}
                dp            = dp + ((3.2780000000000006e-07)*phi*phi +              (6.7084e-07)*phi +  (-4.390000000000004e-05))*pp*pp + ((-7.230999999999999e-06)*phi*phi +            (-2.37482e-05)*phi +  (0.0004909000000000007))*pp +   ((3.285299999999999e-05)*phi*phi +            (9.63723e-05)*phi +               (-0.00115));
            }
            if(sec == 2){
                dp            =               ((-7.741952e-06)*phi*phi + (-2.2402167000000004e-05)*phi + (-0.00042652900000000004))*pp*pp +            ((7.54079e-05)*phi*phi + (-1.3333999999999984e-05)*phi +  (0.0002420100000000004))*pp +            ((-0.000147876)*phi*phi +             (0.00057905)*phi +              (0.0253551));
                if(pp < 7){dp = dp +             ((5.3611e-06)*phi*phi +              (8.1979e-06)*phi +              (5.9789e-04))*pp*pp +            ((-4.8185e-05)*phi*phi +             (-1.5188e-04)*phi +             (-0.0084675))*pp +              ((9.2324e-05)*phi*phi +             (6.4420e-04)*phi +               (0.026792));}
                else{      dp = dp +            ((-6.1139e-05)*phi*phi +              (5.4087e-06)*phi +              (-0.0021284))*pp*pp +              ((0.0010007)*phi*phi +              (9.3492e-05)*phi +               (0.039813))*pp +              ((-0.0040434)*phi*phi +             (-0.0010953)*phi +               (-0.18112));}
                dp            = dp +           ((6.221217e-07)*phi*phi +  (1.9596000000000003e-06)*phi +              (-9.826e-05))*pp*pp +           ((-1.28576e-05)*phi*phi +            (-4.36589e-05)*phi +             (0.00130342))*pp +             ((5.80399e-05)*phi*phi +            (0.000215388)*phi + (-0.0040414000000000005));
            }
            if(sec == 3){
                dp            =      ((-5.115364000000001e-06)*phi*phi + (-1.1983000000000004e-05)*phi +  (-0.0006832899999999999))*pp*pp +            ((4.52287e-05)*phi*phi +  (0.00020855000000000003)*phi +  (0.0034986999999999996))*pp +  ((-9.044610000000001e-05)*phi*phi +            (-0.00106657)*phi +   (0.017954199999999997));
                if(pp < 7){dp = dp +             ((9.9281e-07)*phi*phi +              (3.4879e-06)*phi +               (0.0011673))*pp*pp +            ((-2.0071e-05)*phi*phi +             (-3.1362e-05)*phi +              (-0.012329))*pp +              ((6.9463e-05)*phi*phi +             (3.5102e-05)*phi +               (0.037505));}
                else{      dp = dp +            ((-3.2178e-06)*phi*phi +              (4.0630e-05)*phi +               (-0.005209))*pp*pp +             ((2.0884e-05)*phi*phi +             (-6.8800e-04)*phi +               (0.086513))*pp +              ((3.9530e-05)*phi*phi +              (0.0029306)*phi +                (-0.3507));}
                dp            = dp + ((-4.045999999999999e-07)*phi*phi + (-1.3115999999999994e-06)*phi +  (3.9510000000000006e-05))*pp*pp +              ((5.521e-06)*phi*phi +  (2.4436999999999997e-05)*phi +             (-0.0016887))*pp + ((-1.0962999999999997e-05)*phi*phi +           (-0.000151944)*phi +   (0.009313599999999998));
            }
            if(sec == 4){
                dp            =     ((-3.9278116999999996e-06)*phi*phi +  (2.2289300000000004e-05)*phi +  (0.00012665000000000002))*pp*pp + ((4.8649299999999995e-05)*phi*phi +             (-0.00012554)*phi +  (-0.005955500000000001))*pp + ((-0.00014617199999999997)*phi*phi +            (-0.00028571)*phi +              (0.0606998));
                if(pp < 7){dp = dp +            ((-4.8455e-06)*phi*phi +             (-1.2074e-05)*phi +               (0.0013221))*pp*pp +             ((3.2207e-05)*phi*phi +              (1.3144e-04)*phi +              (-0.010451))*pp +             ((-3.7365e-05)*phi*phi +            (-4.2344e-04)*phi +               (0.019952));}
                else{      dp = dp +            ((-3.9554e-05)*phi*phi +              (5.5496e-06)*phi +              (-0.0058293))*pp*pp +             ((6.5077e-04)*phi*phi +              (2.6735e-05)*phi +               (0.095025))*pp +              ((-0.0026457)*phi*phi +            (-6.1394e-04)*phi +                (-0.3793));}
                dp            = dp +          ((-4.593089e-07)*phi*phi +             (1.40673e-05)*phi +                (6.69e-05))*pp*pp +             ((4.0239e-06)*phi*phi +            (-0.000180863)*phi + (-0.0008272199999999999))*pp + ((-5.1310000000000005e-06)*phi*phi +             (0.00049748)*phi +             (0.00255231));
            }
            if(sec == 5){
                dp            =       ((8.036599999999999e-07)*phi*phi +             (2.58072e-05)*phi +             (0.000360217))*pp*pp + ((-9.932400000000002e-06)*phi*phi +           (-0.0005168531)*phi +              (-0.010904))*pp +  ((1.8516299999999998e-05)*phi*phi +  (0.0015570900000000001)*phi +               (0.066493));
                if(pp < 7){dp = dp +             ((7.7156e-07)*phi*phi +             (-3.9566e-05)*phi +             (-2.3589e-04))*pp*pp +            ((-9.8309e-06)*phi*phi +              (3.7353e-04)*phi +              (0.0020382))*pp +              ((2.9506e-05)*phi*phi +            (-8.0409e-04)*phi +             (-0.0045615));}
                else{      dp = dp +            ((-3.2410e-05)*phi*phi +             (-4.3301e-05)*phi +              (-0.0028742))*pp*pp +             ((5.3787e-04)*phi*phi +              (6.8921e-04)*phi +               (0.049578))*pp +              ((-0.0021955)*phi*phi +             (-0.0027698)*phi +               (-0.21142));}
                dp            = dp +            ((-1.2151e-06)*phi*phi +             (-8.5087e-06)*phi +               (4.968e-05))*pp*pp +            ((1.46998e-05)*phi*phi +             (0.000115047)*phi +            (-0.00039269))*pp + ((-4.0368600000000005e-05)*phi*phi +            (-0.00037078)*phi +             (0.00073998));
            }
            if(sec == 6){
                dp            =     ((-1.9552099999999998e-06)*phi*phi +   (8.042199999999997e-06)*phi + (-2.1324000000000028e-05))*pp*pp + ((1.6969399999999997e-05)*phi*phi +  (-6.306600000000001e-05)*phi +            (-0.00485568))*pp +             ((-2.7723e-05)*phi*phi + (-6.828400000000003e-05)*phi +              (0.0447535));
                if(pp < 7){dp = dp +            ((-8.2535e-07)*phi*phi +              (9.1433e-06)*phi +              (3.5395e-04))*pp*pp +            ((-3.4272e-06)*phi*phi +             (-1.3012e-04)*phi +             (-0.0030724))*pp +              ((4.9211e-05)*phi*phi +             (4.5807e-04)*phi +              (0.0058932));}
                else{      dp = dp +            ((-4.9760e-05)*phi*phi +             (-7.2903e-05)*phi +              (-0.0020453))*pp*pp +             ((8.0918e-04)*phi*phi +               (0.0011688)*phi +               (0.037042))*pp +              ((-0.0032504)*phi*phi +             (-0.0046169)*phi +               (-0.16331));}
                dp            = dp + ((-7.153000000000002e-07)*phi*phi +             (1.62859e-05)*phi +               (8.129e-05))*pp*pp + ((7.2249999999999994e-06)*phi*phi +            (-0.000178946)*phi + (-0.0009485399999999999))*pp + ((-1.3018000000000003e-05)*phi*phi + (0.00046643000000000005)*phi +             (0.00266508));
            }
        }
    }
    //==============================//  Electron Corrections (End)  //==============================//
    
    //==============================//        π+ Corrections        //==============================//
    if(corPip != 0 && ivec == 1){
        if(corPip == 1){ // Fall 2018 - Pass 1 Corrections
            if(sec == 1){
                dp =      ((-5.4904e-07)*phi*phi + (-1.4436e-05)*phi +  (3.1534e-04))*pp*pp +  ((3.8231e-06)*phi*phi +  (3.6582e-04)*phi +  (-0.0046759))*pp + ((-5.4913e-06)*phi*phi + (-4.0157e-04)*phi +    (0.010767));
                dp = dp +  ((6.1103e-07)*phi*phi +  (5.5291e-06)*phi + (-1.9120e-04))*pp*pp + ((-3.2300e-06)*phi*phi +  (1.5377e-05)*phi +  (7.5279e-04))*pp +  ((2.1434e-06)*phi*phi + (-6.9572e-06)*phi + (-7.9333e-05));
                dp = dp + ((-1.3049e-06)*phi*phi +  (1.1295e-05)*phi +  (4.5797e-04))*pp*pp +  ((9.3122e-06)*phi*phi + (-5.1074e-05)*phi +  (-0.0030757))*pp + ((-1.3102e-05)*phi*phi +  (2.2153e-05)*phi +   (0.0040938));
            }
            if(sec == 2){
                dp =      ((-1.0087e-06)*phi*phi +  (2.1319e-05)*phi +  (7.8641e-04))*pp*pp +  ((6.7485e-06)*phi*phi +  (7.3716e-05)*phi +  (-0.0094591))*pp + ((-1.1820e-05)*phi*phi + (-3.8103e-04)*phi +    (0.018936));
                dp = dp +  ((8.8155e-07)*phi*phi + (-2.8257e-06)*phi + (-2.6729e-04))*pp*pp + ((-5.4499e-06)*phi*phi +  (3.8397e-05)*phi +   (0.0015914))*pp +  ((6.8926e-06)*phi*phi + (-5.9386e-05)*phi +  (-0.0021749));
                dp = dp + ((-2.0147e-07)*phi*phi +  (1.1061e-05)*phi +  (3.8827e-04))*pp*pp +  ((4.9294e-07)*phi*phi + (-6.0257e-05)*phi +  (-0.0022087))*pp +  ((9.8548e-07)*phi*phi +  (5.9047e-05)*phi +   (0.0022905));
            }
            if(sec == 3){
                dp =       ((8.6722e-08)*phi*phi + (-1.7975e-05)*phi +  (4.8118e-05))*pp*pp +  ((2.6273e-06)*phi*phi +  (3.1453e-05)*phi +  (-0.0015943))*pp + ((-6.4463e-06)*phi*phi + (-5.8990e-05)*phi +   (0.0041703));
                dp = dp +  ((9.6317e-07)*phi*phi + (-1.7659e-06)*phi + (-8.8318e-05))*pp*pp + ((-5.1346e-06)*phi*phi +  (8.3318e-06)*phi +  (3.7723e-04))*pp +  ((3.9548e-06)*phi*phi + (-6.9614e-05)*phi +  (2.1393e-04));
                dp = dp +  ((5.6438e-07)*phi*phi +  (8.1678e-06)*phi + (-9.4406e-05))*pp*pp + ((-3.9074e-06)*phi*phi + (-6.5174e-05)*phi +  (5.4218e-04))*pp +  ((6.3198e-06)*phi*phi +  (1.0611e-04)*phi + (-4.5749e-04));
            }
            if(sec == 4){
                dp =       ((4.3406e-07)*phi*phi + (-4.9036e-06)*phi +  (2.3064e-04))*pp*pp +  ((1.3624e-06)*phi*phi +  (3.2907e-05)*phi +  (-0.0034872))*pp + ((-5.1017e-06)*phi*phi +  (2.4593e-05)*phi +   (0.0092479));
                dp = dp +  ((6.0218e-07)*phi*phi + (-1.4383e-05)*phi + (-3.1999e-05))*pp*pp + ((-1.1243e-06)*phi*phi +  (9.3884e-05)*phi + (-4.1985e-04))*pp + ((-1.8808e-06)*phi*phi + (-1.2222e-04)*phi +   (0.0014037));
                dp = dp + ((-2.5490e-07)*phi*phi + (-8.5120e-07)*phi +  (7.9109e-05))*pp*pp +  ((2.5879e-06)*phi*phi +  (8.6108e-06)*phi + (-5.1533e-04))*pp + ((-4.4521e-06)*phi*phi + (-1.7012e-05)*phi +  (7.4848e-04));
            }
            if(sec == 5){
                dp =       ((2.4292e-07)*phi*phi +  (8.8741e-06)*phi +  (2.9482e-04))*pp*pp +  ((3.7229e-06)*phi*phi +  (7.3215e-06)*phi +  (-0.0050685))*pp + ((-1.1974e-05)*phi*phi + (-1.3043e-04)*phi +   (0.0078836));
                dp = dp +  ((1.0867e-06)*phi*phi + (-7.7630e-07)*phi + (-4.4930e-05))*pp*pp + ((-5.6564e-06)*phi*phi + (-1.3417e-05)*phi +  (2.5224e-04))*pp +  ((6.8460e-06)*phi*phi +  (9.0495e-05)*phi + (-4.6587e-04));
                dp = dp +  ((8.5720e-07)*phi*phi + (-6.7464e-06)*phi + (-4.0944e-05))*pp*pp + ((-4.7370e-06)*phi*phi +  (5.8808e-05)*phi +  (1.9047e-04))*pp +  ((5.7404e-06)*phi*phi + (-1.1105e-04)*phi + (-1.9392e-04));
            }
            if(sec == 6){
                dp =       ((2.1191e-06)*phi*phi + (-3.3710e-05)*phi +  (2.5741e-04))*pp*pp + ((-1.2915e-05)*phi*phi +  (2.3753e-04)*phi + (-2.6882e-04))*pp +  ((2.2676e-05)*phi*phi + (-2.3115e-04)*phi +   (-0.001283));
                dp = dp +  ((6.0270e-07)*phi*phi + (-6.8200e-06)*phi +  (1.3103e-04))*pp*pp + ((-1.8745e-06)*phi*phi +  (3.8646e-05)*phi + (-8.8056e-04))*pp +  ((2.0885e-06)*phi*phi + (-3.4932e-05)*phi +  (4.5895e-04));
                dp = dp +  ((4.7349e-08)*phi*phi + (-5.7528e-06)*phi + (-3.4097e-06))*pp*pp +  ((1.7731e-06)*phi*phi +  (3.5865e-05)*phi + (-5.7881e-04))*pp + ((-9.7008e-06)*phi*phi + (-4.1836e-05)*phi +   (0.0035403));
            }
        }
        
        if(corPip == 2){ // Spring 2019 - Pass 2 Corrections
            if(sec == 1){
                dp =                   ((1.07338e-06)*phi*phi + (0.00011237500000000001)*phi + (0.00046984999999999996))*pp*pp + ((-2.9323999999999997e-06)*phi*phi + (-0.000777199)*phi + (-0.0061279))*pp + ((3.7362e-06)*phi*phi + (0.00049608)*phi + (0.0156802));
                if(pp < 3.5){dp = dp + ((-8.0699e-06)*phi*phi + (3.3838e-04)*phi + (0.0051143))*pp*pp + ((3.0234e-05)*phi*phi + (-0.0015167)*phi + (-0.023081))*pp + ((-1.3818e-05)*phi*phi + (0.0011894)*phi + (0.015812));}
                else{        dp = dp +  ((2.8904e-07)*phi*phi + (-1.0534e-04)*phi + (-0.0023996))*pp*pp + ((2.3276e-06)*phi*phi + (0.0010502)*phi + (0.022682))*pp + ((-1.9319e-05)*phi*phi + (-0.0025179)*phi + (-0.050285));}
            }
            if(sec == 2){
                dp =                   ((2.97335e-06)*phi*phi + (7.68257e-05)*phi + (0.001132483))*pp*pp + ((-1.86553e-05)*phi*phi + (-0.000511963)*phi + (-0.0111051))*pp + ((2.16081e-05)*phi*phi + (0.000100984)*phi + (0.0189673));
                if(pp < 3.5){dp = dp + ((-1.4761e-06)*phi*phi + (4.9397e-06)*phi + (0.0014986))*pp*pp + ((6.4311e-06)*phi*phi + (-3.8570e-05)*phi + (-0.005309))*pp + ((2.2896e-06)*phi*phi + (-1.8426e-04)*phi + (-0.0030622));}
                else{        dp = dp +  ((3.3302e-06)*phi*phi + (-8.4794e-05)*phi + (-0.0020262))*pp*pp + ((-3.5962e-05)*phi*phi + (9.1367e-04)*phi + (0.019333))*pp + ((9.5116e-05)*phi*phi + (-0.0023371)*phi + (-0.045778));}
            }
            if(sec == 3){
                dp =        ((1.9689700000000002e-07)*phi*phi + (-6.73721e-05)*phi + (0.001145664))*pp*pp + ((-1.3357999999999998e-07)*phi*phi + (0.0004974620000000001)*phi + (-0.01087555))*pp + ((5.23389e-06)*phi*phi + (-0.00038631399999999996)*phi + (0.012021909999999999));
                if(pp < 3.5){dp = dp + ((-3.7071e-06)*phi*phi + (-6.7985e-05)*phi + (0.0073195))*pp*pp + ((1.2081e-05)*phi*phi + (4.0719e-04)*phi + (-0.032716))*pp + ((1.8109e-06)*phi*phi + (-5.6304e-04)*phi + (0.022124));}
                else{        dp = dp +  ((2.9228e-06)*phi*phi + (-7.4216e-07)*phi + (-0.0033922))*pp*pp + ((-2.7026e-05)*phi*phi + (-7.5709e-06)*phi + (0.03267))*pp + ((5.8592e-05)*phi*phi + (3.8319e-05)*phi + (-0.076661));}
            }
            if(sec == 4){
                dp =                    ((5.4899e-07)*phi*phi + (-1.82236e-05)*phi + (0.0007486388))*pp*pp + ((-1.0743e-06)*phi*phi + (0.000125103)*phi + (-0.00743795))*pp + ((1.9187e-06)*phi*phi + (-5.0545e-05)*phi + (0.01528271));
                if(pp < 3.5){dp = dp + ((-7.1834e-06)*phi*phi + (1.2815e-04)*phi + (0.004323))*pp*pp + ((2.7688e-05)*phi*phi + (-4.9122e-04)*phi + (-0.020112))*pp + ((-1.5879e-05)*phi*phi + (3.5148e-04)*phi + (0.013367));}
                else{        dp = dp + ((-2.2635e-06)*phi*phi + (3.3612e-05)*phi + (-0.0024779))*pp*pp + ((2.7765e-05)*phi*phi + (-4.4868e-04)*phi + (0.02433))*pp + ((-7.6567e-05)*phi*phi + (0.0013553)*phi + (-0.058136));}
            }
            if(sec == 5){
                dp =                    ((9.5628e-07)*phi*phi + (-1.4e-06)*phi + (0.00116279))*pp*pp + ((-3.723047e-06)*phi*phi + (2.09447e-05)*phi + (-0.0101853))*pp + ((9.326299999999999e-06)*phi*phi + (-0.0001111214)*phi + (0.0130134));
                if(pp < 3.5){dp = dp + ((-8.2807e-06)*phi*phi + (-1.2620e-04)*phi + (0.0060821))*pp*pp + ((3.8915e-05)*phi*phi + (6.3989e-04)*phi + (-0.028784))*pp + ((-3.7765e-05)*phi*phi + (-7.0844e-04)*phi + (0.021177));}
                else{        dp = dp + ((-8.7415e-08)*phi*phi + (3.5806e-05)*phi + (-0.0022065))*pp*pp + ((5.3612e-06)*phi*phi + (-4.2740e-04)*phi + (0.022369))*pp + ((-2.3587e-05)*phi*phi + (0.0011096)*phi + (-0.056773));}
            }
            if(sec == 6){
                dp =                   ((5.86478e-07)*phi*phi + (3.5833999999999994e-06)*phi + (0.00108574))*pp*pp + ((-4.433118e-06)*phi*phi + (-5.3565999999999995e-05)*phi + (-0.00873827))*pp + ((2.0270600000000002e-05)*phi*phi + (-7.0902e-05)*phi + (0.0077521));
                if(pp < 3.5){dp = dp +  ((1.4952e-06)*phi*phi + (1.3858e-05)*phi + (0.0028677))*pp*pp + ((-8.0852e-06)*phi*phi + (-1.1384e-04)*phi + (-0.015643))*pp + ((9.5078e-06)*phi*phi + (1.3285e-04)*phi + (0.014019));}
                else{        dp = dp + ((-5.7308e-07)*phi*phi + (-3.8697e-05)*phi + (-0.0030495))*pp*pp + ((1.0905e-05)*phi*phi + (3.8288e-04)*phi + (0.030355))*pp + ((-3.1873e-05)*phi*phi + (-9.6019e-04)*phi + (-0.074345));}
            }
        }
        
        if(corPip == 3){ // Fall 2018 - Pass 2 Corrections
            if(sec == 1){
                dp              =           ((1.338454e-06)*phi*phi +   (4.714629999999999e-05)*phi +  (0.00014719))*pp*pp + ((-2.8460000000000004e-06)*phi*phi +            (-0.000406925)*phi +           (-0.00367325))*pp +           ((-1.193548e-05)*phi*phi +            (-0.000225083)*phi +           (0.01544091));
                if(pp < 2.5){dp = dp +        ((1.0929e-05)*phi*phi +             (-3.8002e-04)*phi +    (-0.01412))*pp*pp +             ((-2.8491e-05)*phi*phi +              (5.0952e-04)*phi +              (0.037728))*pp +              ((1.6927e-05)*phi*phi +              (1.8165e-04)*phi +            (-0.027772));}
                else{        dp = dp +        ((4.3191e-07)*phi*phi +             (-9.0581e-05)*phi +  (-0.0011766))*pp*pp +             ((-3.6232e-06)*phi*phi +               (0.0010342)*phi +              (0.012454))*pp +              ((1.2235e-05)*phi*phi +              (-0.0025855)*phi +            (-0.035323));}
                dp              = dp +       ((-3.7494e-07)*phi*phi +             (-1.5439e-06)*phi +  (4.2760e-05))*pp*pp +              ((3.5348e-06)*phi*phi +              (4.8165e-05)*phi +           (-2.3799e-04))*pp +             ((-8.2116e-06)*phi*phi +             (-7.1750e-05)*phi +           (1.5984e-04));
            }
            if(sec == 2){
                dp              =             ((5.8222e-07)*phi*phi +  (5.0666599999999994e-05)*phi +  (0.00051782))*pp*pp +              ((3.3785e-06)*phi*phi +            (-0.000343093)*phi + (-0.007453400000000001))*pp + ((-2.2014899999999998e-05)*phi*phi + (-0.00027579899999999997)*phi + (0.015119099999999998));
                if(pp < 2.5){dp = dp +        ((9.2373e-06)*phi*phi +             (-3.3151e-04)*phi +   (-0.019254))*pp*pp +             ((-2.7546e-05)*phi*phi +              (5.3915e-04)*phi +              (0.052516))*pp +              ((2.5220e-05)*phi*phi +              (7.5362e-05)*phi +            (-0.033504));}
                else{        dp = dp +        ((2.2654e-08)*phi*phi +             (-8.8436e-05)*phi +  (-0.0013542))*pp*pp +              ((3.0630e-07)*phi*phi +              (9.4319e-04)*phi +                (0.0147))*pp +             ((-3.5941e-06)*phi*phi +              (-0.0022473)*phi +            (-0.036874));}
                dp              = dp +        ((4.3694e-07)*phi*phi +              (1.1476e-05)*phi +  (1.1123e-04))*pp*pp +             ((-2.4617e-06)*phi*phi +             (-7.5353e-05)*phi +           (-6.2511e-04))*pp +             ((-1.0387e-06)*phi*phi +              (5.8447e-05)*phi +           (6.4986e-04));
            }
            if(sec == 3){
                dp              =           ((-6.17815e-07)*phi*phi + (-1.4503600000000001e-05)*phi + (0.000584689))*pp*pp +             ((8.27871e-06)*phi*phi +              (9.2796e-05)*phi +         (-0.0078185692))*pp + ((-1.6866360000000002e-05)*phi*phi +  (-8.065000000000001e-05)*phi +            (0.0159476));
                if(pp < 2.5){dp = dp +        ((1.8595e-06)*phi*phi +              (3.6900e-04)*phi +  (-0.0099622))*pp*pp +              ((8.4410e-06)*phi*phi +              (-0.0010457)*phi +              (0.027038))*pp +             ((-1.2191e-05)*phi*phi +              (6.0203e-04)*phi +            (-0.019176));}
                else{        dp = dp +        ((6.8265e-07)*phi*phi +              (3.0246e-05)*phi +  (-0.0011116))*pp*pp +             ((-4.8481e-06)*phi*phi +             (-3.7082e-04)*phi +              (0.011452))*pp +              ((7.2478e-06)*phi*phi +              (9.9858e-04)*phi +            (-0.027972));}
                dp              = dp +        ((1.8639e-07)*phi*phi +              (4.9444e-06)*phi + (-2.9030e-05))*pp*pp +             ((-1.3752e-06)*phi*phi +             (-3.3709e-05)*phi +            (3.8288e-04))*pp +              ((1.0113e-06)*phi*phi +              (5.1273e-05)*phi +          (-6.7844e-04));
            }
            if(sec == 4){
                dp              =  ((9.379499999999998e-07)*phi*phi + (-2.8101700000000002e-05)*phi +  (0.00053373))*pp*pp + ((-1.6185199999999991e-06)*phi*phi +  (0.00017444500000000001)*phi + (-0.005648269999999999))*pp +  ((-3.495700000000003e-06)*phi*phi +  (-7.845739999999999e-05)*phi + (0.010768400000000001));
                if(pp < 2.5){dp = dp +        ((9.5779e-06)*phi*phi +              (3.5339e-04)*phi +    (-0.01054))*pp*pp +             ((-1.8077e-05)*phi*phi +              (-0.0010543)*phi +              (0.028379))*pp +              ((3.1773e-06)*phi*phi +              (5.6223e-04)*phi +            (-0.018865));}
                else{        dp = dp +        ((7.7000e-07)*phi*phi +              (4.1000e-06)*phi +  (-0.0010144))*pp*pp +             ((-8.1960e-06)*phi*phi +             (-4.7753e-05)*phi +              (0.010594))*pp +              ((2.0716e-05)*phi*phi +              (1.2151e-04)*phi +            (-0.028619));}
                dp              = dp +        ((4.8394e-07)*phi*phi +              (3.6342e-06)*phi + (-2.0136e-04))*pp*pp +             ((-3.2757e-06)*phi*phi +             (-3.5397e-05)*phi +             (0.0015599))*pp +              ((3.2095e-06)*phi*phi +              (7.9013e-05)*phi +            (-0.002012));
            }
            if(sec == 5){
                dp              = ((1.7566900000000006e-07)*phi*phi +             (2.21337e-05)*phi +   (0.0011632))*pp*pp +   ((2.812770000000001e-06)*phi*phi + (-0.00018654499999999998)*phi + (-0.011854620000000001))*pp +  ((-8.442900000000003e-06)*phi*phi + (-0.00011505800000000001)*phi +            (0.0176174));
                if(pp < 2.5){dp = dp +        ((3.3685e-05)*phi*phi +              (2.8972e-04)*phi +   (-0.017862))*pp*pp +             ((-8.4089e-05)*phi*phi +             (-9.8038e-04)*phi +              (0.050405))*pp +              ((4.3478e-05)*phi*phi +              (6.9924e-04)*phi +            (-0.033066));}
                else{        dp = dp +        ((4.6106e-07)*phi*phi +             (-3.6786e-05)*phi +  (-0.0015894))*pp*pp +             ((-4.4217e-06)*phi*phi +              (3.7321e-04)*phi +              (0.015917))*pp +              ((7.5188e-06)*phi*phi +             (-8.0676e-04)*phi +            (-0.036944));}
                dp              = dp +        ((4.3113e-07)*phi*phi +              (2.6869e-06)*phi + (-2.1326e-04))*pp*pp +             ((-3.1063e-06)*phi*phi +             (-2.7152e-05)*phi +             (0.0017964))*pp +              ((3.1946e-06)*phi*phi +              (4.2059e-05)*phi +           (-0.0031325));
            }
            if(sec == 6){
                dp              =            ((1.94354e-06)*phi*phi +  (1.3306000000000006e-05)*phi +  (0.00067634))*pp*pp +             ((-7.9584e-06)*phi*phi +  (-7.949999999999998e-05)*phi + (-0.005861990000000001))*pp +   ((6.994000000000005e-07)*phi*phi +             (-0.00022435)*phi +            (0.0118564));
                if(pp < 2.5){dp = dp +        ((1.7381e-05)*phi*phi +              (5.4630e-04)*phi +   (-0.019637))*pp*pp +             ((-3.8681e-05)*phi*phi +              (-0.0017358)*phi +                (0.0565))*pp +              ((1.2268e-05)*phi*phi +               (0.0011412)*phi +            (-0.035608));}
                else{        dp = dp +       ((-8.9398e-08)*phi*phi +             (-1.2347e-05)*phi +  (-0.0018442))*pp*pp +              ((7.8164e-08)*phi*phi +              (1.3063e-04)*phi +               (0.01783))*pp +              ((8.2374e-06)*phi*phi +             (-3.5862e-04)*phi +            (-0.047011));}
                dp              = dp +        ((4.9123e-07)*phi*phi +              (5.1828e-06)*phi + (-1.3898e-04))*pp*pp +             ((-3.4108e-06)*phi*phi +             (-5.0009e-05)*phi +             (0.0014879))*pp +              ((4.0320e-06)*phi*phi +              (6.5853e-05)*phi +           (-0.0032227));
            }
            
        }
    }
    //==============================//     π+ Corrections (End)     //==============================//
    
    //==============================//        π- Corrections        //==============================//
    if(corPim != 0 && ivec == 2){
        if(sec == 1){ // Fall 2018 - Pass 1 Corrections (Only)
            dp =      ((-4.0192658422317425e-06)*phi*phi -  (2.660222128967742e-05)*phi + 0.004774434682983547)*pp*pp;
            dp = dp +  ((1.9549520962477972e-05)*phi*phi -    0.0002456062756770577*phi - 0.03787692408323466)*pp; 
            dp = dp +   (-2.128953094937459e-05)*phi*phi +    0.0002461708852239913*phi + 0.08060704449822174 - 0.01;
        }
        if(sec == 2){
            dp =        ((1.193010521758372e-05)*phi*phi -  (5.996221756031922e-05)*phi + 0.0009093437955814359)*pp*pp;
            dp = dp +   ((-4.89113824430594e-05)*phi*phi +   0.00021676479488147118*phi - 0.01861892053916726)*pp;  
            dp = dp +    (4.446394152208071e-05)*phi*phi - (3.6592784167335244e-05)*phi + 0.05498710249944096 - 0.01;
        }
        if(sec == 3){
            dp =      ((-1.6596664895992133e-07)*phi*phi +  (6.317189710683516e-05)*phi + 0.0016364212312654086)*pp*pp;
            dp = dp +  ((-2.898409777520318e-07)*phi*phi -   0.00014531513577533802*phi - 0.025456145839203827)*pp;  
            dp = dp +   (2.6432552410603506e-06)*phi*phi +   0.00018447151306275443*phi + 0.06442602664627255 - 0.01;
        }
        if(sec == 4){
            dp =       ((2.4035259647558634e-07)*phi*phi -  (8.649647351491232e-06)*phi + 0.004558993439848128)*pp*pp;
            dp = dp +  ((-5.981498144060984e-06)*phi*phi +   0.00010582131454222416*phi - 0.033572004651981686)*pp;  
            dp = dp +     (8.70140266889548e-06)*phi*phi -   0.00020137414379966883*phi + 0.07258774523336173 - 0.01;   
        }
        if(sec == 5){
            dp =       ((2.5817024702834863e-06)*phi*phi +   0.00010132810066914441*phi + 0.003397314538804711)*pp*pp;
            dp = dp + ((-1.5116941263931812e-05)*phi*phi -   0.00040679799541839254*phi - 0.028144285760769876)*pp;  
            dp = dp +   (1.4701931057951464e-05)*phi*phi +    0.0002426350390593454*phi + 0.06781682510174941 - 0.01;
        }
        if(sec == 6){
            dp =       ((-8.196823669099362e-07)*phi*phi -  (5.280412421933636e-05)*phi + 0.0018457238328451137)*pp*pp;
            dp = dp +  ((5.2675062282094536e-06)*phi*phi +    0.0001515803461044587*phi - 0.02294371578470564)*pp;  
            dp = dp +   (-9.459454671739747e-06)*phi*phi -    0.0002389523716779765*phi + 0.06428970810739926 - 0.01;
        }
    }
    //==============================//     π- Corrections (End)     //==============================//
    
    //==============================//      Proton Corrections      //==============================//
    if(corPro != 0 && ivec == 3){
        if(sec == 1){ // Fall 2018 - Pass 1 Corrections (Only)
            dp = ((1 + TMath::Sign(1, (pp - 1.4)))/2)*((4.4034e-03)*pp   + (-0.01703))    + ((1 + TMath::Sign(1, -(pp - 1.4)))/2)*((-0.10898)*(pp  - 1.4)*(pp  - 1.4)  + (-0.09574)*(pp - 1.4)  + ((4.4034e-03)*1.4   + (-0.01703)));
        }
        if(sec == 2){
            dp = ((1 + TMath::Sign(1, (pp - 1.5)))/2)*((0.01318)*pp      + (-0.03403))    + ((1 + TMath::Sign(1, -(pp - 1.5)))/2)*((-0.09829)*(pp  - 1.5)*(pp  - 1.5)  +  (-0.0986)*(pp - 1.5)  + ((0.01318)*1.5      + (-0.03403)));
        }
        if(sec == 3){
            dp = ((1 + TMath::Sign(1, (pp - 1.05)))/2)*((-4.7052e-03)*pp + (1.2410e-03))  + ((1 + TMath::Sign(1, -(pp - 1.05)))/2)*((-0.22721)*(pp - 1.05)*(pp - 1.05) + (-0.09702)*(pp - 1.05) + ((-4.7052e-03)*1.05 + (1.2410e-03)));
        }
        if(sec == 4){
            dp = ((1 + TMath::Sign(1, (pp - 1.4)))/2)*((-1.0900e-03)*pp  + (-4.0573e-03)) + ((1 + TMath::Sign(1, -(pp - 1.4)))/2)*((-0.09236)*(pp  - 1.4)*(pp  - 1.4)  +   (-0.073)*(pp - 1.4)  + ((-1.0900e-03)*1.4  + (-4.0573e-03)));
        }
        if(sec == 5){
            dp = ((1 + TMath::Sign(1, (pp - 1.5)))/2)*((7.3965e-03)*pp   + (-0.02428))    + ((1 + TMath::Sign(1, -(pp - 1.5)))/2)*((-0.09539)*(pp  - 1.5)*(pp  - 1.5)  + (-0.09263)*(pp - 1.5)  + ((7.3965e-03)*1.5   + (-0.02428)));
        }
        if(sec == 6){
            dp = ((1 + TMath::Sign(1, (pp - 1.15)))/2)*((-7.6214e-03)*pp + (8.1014e-03))  + ((1 + TMath::Sign(1, -(pp - 1.15)))/2)*((-0.12718)*(pp - 1.15)*(pp - 1.15) + (-0.06626)*(pp - 1.15) + ((-7.6214e-03)*1.15 + (8.1014e-03)));
        }
    }
    //==============================//   End of Proton Correction   //==============================//

    return dp/pp;
};


//smearing

TRandom3 *myMC;  // Needed for smearing
void fastMC(void){
        myMC = new TRandom3(0);
}


//New Smearing:
auto smear_oneFactor = [](const TLorentzVector V4rec){
	// smear factor:
	
	double elTh_indeg = V4rec.Theta() * 57.2958;
	double sigmaRes = 0;
 	const double fp_4[5] = {     2.3762794461143009e+000,
    -2.8635604070815707e-001,
     3.8197546356590874e-002,
    -1.8427728047602948e-003,
     2.7446233737351622e-005};
     
	const double factor = fp_4[0] + fp_4[1]*elTh_indeg + fp_4[2]*elTh_indeg*elTh_indeg + fp_4[3]*elTh_indeg*elTh_indeg*elTh_indeg + fp_4[4]*elTh_indeg*elTh_indeg*elTh_indeg*elTh_indeg;;
	
	sigmaRes = 0.008007538378779988 -0.0010764451335380787* elTh_indeg + 6.943332514486485e-05* elTh_indeg * elTh_indeg  -1.1181109847470665e-06* elTh_indeg * elTh_indeg * elTh_indeg ;

    // True generated values (i.e., values of the unsmeared TLorentzVector)
    double inM = V4rec.M();
    double gausValue = myMC->Gaus(0,1);
    double smeared_P  = V4rec.P() +  V4rec.P() * sigmaRes * factor * gausValue;
    
//cout<<inM<<endl;
//    cout<<V4rec.P() * sigmaRes * factor * gausValue << " gaus:"<< gausValue <<endl;   
    double smeared_Th = V4rec.Theta();
    double smeared_Phi = V4rec.Phi();
    
    TLorentzVector V4_new(V4rec.X(), V4rec.Y(), V4rec.Z(), V4rec.E());
    
    V4_new.SetE( TMath::Sqrt( smeared_P*smeared_P + inM*inM )  );
    V4_new.SetRho( smeared_P );
    V4_new.SetTheta( smeared_Th );
    V4_new.SetPhi( smeared_Phi );
    
    return V4_new;
};


// sector should be 1:6
inline auto Get4mom_corr(double ex, double ey, double ez, int sec_mom_corr, int isData) {
    double p_squared = ex * ex + ey * ey + ez * ez;

    // Handle zero or near-zero momentum safely
    if (p_squared < 1e-6) {
        return TLorentzVector(0, 0, 0, 0);  // Return a default zero 4-vector
    }

    if (isData == 0) {
        TLorentzVector V4(ex, ey, ez, sqrt(p_squared + m_e * m_e));
        return smear_oneFactor(V4);  // safe since V4 has valid magnitude
    } else {
        auto fe = dppC(ex, ey, ez, (int)lrint(sec_mom_corr), 0, 1, 0, 0, 0) + 1;
        double corrected_p_squared = fe * fe * p_squared;
        double energy = sqrt(corrected_p_squared + m_e * m_e);
        TLorentzVector elec_corrected(fe * ex, fe * ey, fe * ez, energy);
        return elec_corrected;
    }
}


struct line{
    double k;
    double b;
};
double func(double x, double k, double b){
    return k * x + b;
}
double isOutOfLines(double x, double y, line topLine, line botLine){
    return y > func(x, topLine.k, topLine.b) || y < func(x, botLine.k, botLine.b);
}

double isBetweenOfLines(double x, double y, line topLine, line botLine){
    return y < func(x, topLine.k, topLine.b) && y > func(x, botLine.k, botLine.b);
}

//cutLevel = 1 should be used in all the cuts. O,2 were used for systematic studies to vary the cuts.

//Sectors# are 0-5 if I do not spicify otherwise


// theta and phi in deg:
//float toRD = 57.2958; // 360/2pi
//float thetaEMeasured = p4_ele[i].Theta()*toRD;

// To get Phi you do:
//float phiDC = p4_ele[i].Phi()*toRD;
//if (phiDC < 0 && sectorElectron > 0)  phiDC = phiDC + 360. - sectorElectron * 60.;
// else phiDC=phiDC - sectorElectron * 60.;

// Cut out an odd structure in pass-1 should not be needed for pass-2



bool phiSpikeCut(const double phi, const double theta, int sec, const int cutLevel){
    sec = sec-1;
	double shiftSys = 0;
	if (cutLevel == 0) shiftSys = 1.;
	if (cutLevel == 2) shiftSys = -1.;
    
    // Sec 1:
    if (sec == 0){
        if (theta > (13.5 + shiftSys) && theta < (27 - shiftSys))
            if (phi > (-5 + shiftSys) && phi < (0 - shiftSys))
                return false;
    }
    // Sec 4:
    if (sec == 3){
        if (theta > (13.5 + shiftSys) && theta < (29.5 - shiftSys))
            if (phi > (-6 + shiftSys)&& phi < (-1 - shiftSys))
                return false;
    }
    return true; 
    
}


// Cut on Vertex, 
//cutLevel = 1 should be used
//use   out_tree.Branch("p4_ele_vz", &p4_ele_vz);
bool CutVz(const float _vz, const int cutLevel){
    float minVz[3] = {-8.5, -8, -7.5};
    float maxVz[3] = {2.5,2,1.5};

    return (_vz > minVz[cutLevel] && _vz < maxVz[cutLevel]);
}

// Cut to remove bad PMTs from CALs
// SECTOR IS 1:6, DO NOT FORGRT TO ADD TO SECTOR VARIABLE!!!!!
// Cut level = 1
//  out_tree.Branch("pcalHX", &pcalHX);
//  out_tree.Branch("pcalHY", &pcalHY);
//  out_tree.Branch("pcalHZ", &pcalHZ);
// and etc.

bool BadElementKnockOut(double Hx_pcal, double Hy_pcal, double Hx_ecin, double Hy_ecin, double Hx_ecout, double Hy_ecout, int sector, int cutLevel){
    
    double widthChange = 0.25;
    if (cutLevel == 0)  widthChange = -0.25;
    if (cutLevel == 2)  widthChange = 0.75;
    
    
    //MANUAL LINES NOT FROM CCDB, I LOOKED DATA/SIM Plots and found abnormal behav.
    
    float pcalManualParams[6][2] = {{107.2766, -10602.9779}, //S2
                                    {98.9667, -10262.0167}, //S2
                                    
        
                                    {98.0644, 5825.4023}, //S5
                                    {99.9337, 5098.3456}, //S5
                                    
                                     //S6
                                    {0.4547, -275.9317},
                                    {0.4547, -285.9317}}; //S6
    
    if (sector == 6){
        double k = tan(-30.6*3.1415/180);
        double b = -185;
        // wire swap remove
        if ((isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange} , {k, b - widthChange - 2}) &&
               isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 8.3} , {k, b - widthChange - 8.3 - 2.2})
               ) == false ) return false;
        
        
        float k1 = pcalManualParams[4][0];
        float b1 = pcalManualParams[4][1];
        float k2 = pcalManualParams[5][0];
        float b2 = pcalManualParams[5][1];
        
        // manual line
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        return true;
    }
    
    
    //PCAL
    //Sec 2, Lay=1, Comp = 24
    //Sec 2, Lay=2,Comp = 43
    //Sec 3, Lay=2,Comp = 6
    //Sec 3, Lay=3, Comp = 41
    //Sec 4, Lay=1, Comp = 60
    //Sec 4, Lay=2, Comp = 20

    float linesParams[8][2] = {{-0.5774, 156.9137}, // REMOVED MAY 2023 Sec 2, Lay=1, Comp = 24
                                {-0.5878, 152.4101}, // REMOVED MAY 2023 Sec 2,
                                {0.5897, 120.7937}, // KEEP MAY 2023 
                                {0.5913, 114.3872}, // KEEP MAY 2023
                                
                                //{-82.3095, -24692.8538}, // -302.38 - -313.71
                                //{114.2353, 36036.7541}, //
                                {-0.5934, 129.9189},// REMOVED MAY 2023 Sec 3, Lay=3, Comp = 41
                                {-0.5899, 125.2137},// REMOVED MAY 2023
                                
                                //{1484.1452, 455714.6153}, // -306.95 - 315.36
                                //{-2214.6082, -698249.0095}, // 
                                {0.7971, 288.6657},// REMOVED MAY 2023 Sec 4, Lay=2, Comp = 20
                                {0.7848, 277.9871}};// REMOVED MAY 2023 Sec 4, Lay=2, Comp = 20
    

    
    //ECIN
    // Sec 1, Lay=5,Comp = 8
    // Sec 1, Lay=6,Comp = 19
    // Sec 1, Lay=6,Comp = 21
    // Sec 4, Lay=5,Comp = 13
    // Sec 4, Lay=5,Comp = 15
    // Sec 5, Lay=5,Comp = 16
    

    // ECOUT
    // Sec 1, Lay=7,Comp = 32
    // Sec 1, Lay=9,Comp = 19
    // Sec 1, Lay=9,Comp = 21
    // Sec 1, Lay=9,Comp = 24
    // Sec 5, Lay=7,Comp = 18
    
    float ecoutParams[8][2] = {//{351.542, -129188.6636}, // replace by 357.14 - 367.99
                                //{-1656.4878, 591767.5446}, // reverse
                                {0.5765, -222.5612},
                                {0.5618, -233.3443},
                                {0.572, -199.4151},
                                {0.5688, -210.0406},
                                {0.5684, -165.1428},
                                {0.5669, -175.8632},
                                {-0.5841, -252.1105}, // KEEP
                                {-0.5775, -263.2072}}; // KEEP
                                
    
    // SECTOR 1
    if (sector == 1){
        
        double k = tan(29.5*3.1415/180);
        double b = -92;
        if ( (isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange} , {k, b - widthChange - 2.4}) &&
                isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 9.1} , {k, b - widthChange - 9.1 - 2.4}) &&
                isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 127} , {k, b - widthChange - 127 - 2.4}) &&
                isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange - 127 - 8} , {k, b - widthChange -127 - 8 - 2.4}) 
               ) == false) return false;
        
        
        return true;
        
    }

    
    
    // SECTOR 2
    if (sector == 2){
     
        float k1 = linesParams[2][0];
        float b1 = linesParams[2][1];
        float k2 = linesParams[3][0];
        float b2 = linesParams[3][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        k1 = pcalManualParams[0][0];
        b1 = pcalManualParams[0][1];
        k2 = pcalManualParams[1][0];
        b2 = pcalManualParams[1][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        return true;
        
    }
    
    // SECTOR 3
    if (sector == 3){
    
    	if (Hx_pcal > -313.71 && Hx_pcal <  -302.38) return false;
        
        return true;
        
    }
    
    // SECTOR 4
    if (sector == 4){
    
        double k = tan(-29.6*3.1415/180);
        double b = -232.8;
        
        if (Hx_pcal > -127.5 && Hx_pcal < -122.5) return false;
        return (isOutOfLines(Hx_pcal, Hy_pcal, {k, b + widthChange} , {k, b - widthChange - 3.5}));
    
        
        return true;
        
    }
    
    // SECTOR 5
    if (sector == 5){
        
        float k1 = ecoutParams[6][0];
        float b1 = ecoutParams[6][1];
        float k2 = ecoutParams[7][0];
        float b2 = ecoutParams[7][1];
        
        //ecout July 2023
        if (isBetweenOfLines(Hx_ecout, Hy_ecout,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
       
        k1 = pcalManualParams[2][0];
        b1 = pcalManualParams[2][1];
        k2 = pcalManualParams[3][0];
        b2 = pcalManualParams[3][1];
        
        if (isBetweenOfLines(Hx_pcal, Hy_pcal,  {k1, b1 + widthChange}, {k2, b2 - widthChange})) return false;
        
        return true;
        
    }
    
    return true;
}


TVector3 rotVect(TVector3 resVector3, const int sec){
	resVector3.RotateZ(-60*sec/57.2958);
	resVector3.RotateY(-25/57.2958);
    return resVector3;
}


//Var. come from:
//  out_tree.Branch("dcXR1", &dcXR1);
//  out_tree.Branch("dcYR1", &dcYR1);
//  out_tree.Branch("dcYZ1", &dcZR1);
    
//  out_tree.Branch("dcXR2", &dcXR2);
//  out_tree.Branch("dcYR2", &dcYR2);
//  out_tree.Branch("dcZR2", &dcZR2);
    
//  out_tree.Branch("dcXR3", &dcXR3);
//  out_tree.Branch("dcYR3", &dcYR3);
//  out_tree.Branch("dcZR3", &dcZR3);


// This is how you do rotation for DC:

// DC rotating                        
//TVector3 dcR1Vector3(dc_XR1[i], dc_YR1[i], dc_ZR1[i]);
//TVector3 dcR2Vector3(dc_XR2[i], dc_YR2[i], dc_ZR2[i]);
//TVector3 dcR3Vector3(dc_XR3[i], dc_YR3[i], dc_ZR3[i]);
                        
//rotVect(dcR1Vector3, sectorElectron);
//rotVect(dcR2Vector3, sectorElectron);
//rotVect(dcR3Vector3, sectorElectron);
                        
//DCXY dcParam = {dcR1Vector3.X(), dcR1Vector3.Y(),
//                dcR2Vector3.X(), dcR2Vector3.Y(),
//                dcR3Vector3.X(), dcR3Vector3.Y() };
                                        
                                        

///END INPUT FOR DC:
/////////////////////////////////////////////////////////////////////////////////////////////////



bool DCFidXY(float X, float Y, int region, int sector, int cutLevel){
    
    float cutLimit = 0;
    float angle = 0.495;
    
    if (region == 1){
        cutLimit = 72;
        angle = 0.50;
       }
    if (region == 2){
        cutLimit = 114;
        angle = 0.505;
       }
    if (region == 3){
        cutLimit = 180;
	}
	
	cutLimit -= 2*(cutLevel - 1) * region * 0.6;
	
    return (Y >= -angle*(X + cutLimit) && Y <= angle*(X + cutLimit));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// DC FID CUT 
struct DCXY{
 double r1X,r1Y,r2X,r2Y,r3X,r3Y; // rotatesd x and y coordinates
};
bool CutDCfid(const DCXY& dc, int sec, const int cutLevel){
    return (
         DCFidXY(dc.r1X, dc.r1Y, 1, sec + 1, cutLevel) && 
         DCFidXY(dc.r2X, dc.r2Y, 2, sec + 1, cutLevel) &&
         DCFidXY(dc.r3X, dc.r3Y, 3, sec + 1, cutLevel)
         );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///INPUT FOR DC:

////////////////////////////////
// sampling fraction cut:

/// Input for SF cut (deposited energy divided by momentum):
/// sf = kin.calEnerg / p4_electron.P(), Edep = kin.calEnerg, for cal energy, sum all 3
// Edep = kin.calEnerg = ecout_E[i] + ecin_E[i] + pcal_E[i]
// p4_electron is scattered electron.

// Energy comes from:
//  out_tree.Branch("pcalE", &pcalE);
//  out_tree.Branch("ecinE", &ecinE);
//  out_tree.Branch("ecoutE", &ecoutE);

bool SfCutValerii_Edepos(const double sf,const double Edep,const int sec,const int cutLevel,const int isData){
    double SigmaRange = 3.5;
    double sigmaRangeChange = 0.5;
    
    if (cutLevel == 0) SigmaRange += sigmaRangeChange;
        
    if (cutLevel == 2) SigmaRange -= sigmaRangeChange;    
    
    double meanAll[3][7] = {{0.28617, 0.27974, 0.27476, 0.27264, 0.27074, 0.27575, 0.29045},
    						{-0.04047, -0.03786, -0.03409, -0.03304, -0.03211, -0.03394, -0.03974},
    						{-0.00296, -0.00118, -0.00141, -0.00067, 0.00051, -0.00143, -0.00286}};
    						
    double sigmaAll[3][7] = {{0.0172, 0.01878, 0.017, 0.01566, 0.01586, 0.01684, 0.01485},
    						 {-0.00122, -0.00296, -0.00224, 0.00028, -0.00135, -0.00185, -0.00053},
    						 {-0.0012, -0.00135, -0.00129, -0.00131, -0.00092, -0.00116, -0.00137 }};
    
    
    if (isData){

        
        double mean = meanAll[0][sec] + meanAll[1][sec] / Edep + meanAll[2][sec]/ ( Edep * Edep);
        double sigma = sigmaAll[0][sec] + sigmaAll[1][sec] / Edep + sigmaAll[2][sec] / ( Edep * Edep);
        
        double lowSFcut = mean - sigma * SigmaRange;
        
        return sf > lowSFcut;
    }
    else{

        double mean = meanAll[0][6]  + meanAll[1][6]  / Edep + meanAll[2][6] / ( Edep * Edep);
        double sigma = sigmaAll[0][6] + sigmaAll[1][6]  / Edep + sigmaAll[2][6]  / ( Edep * Edep);
        
        double lowSFcut = mean - sigma * SigmaRange;
        
        return sf > lowSFcut;
    }
    return 0;
}


// PCAL FID CUT 
// uses 
//  out_tree.Branch("pcalLu", &pcalLu);
//  out_tree.Branch("pcalLv", &pcalLv);
//  out_tree.Branch("pcalLw", &pcalLw);    

		bool PCALFid_VW(float V, float W, float U, int cutLevel){
		
			if (cutLevel < 0 || cutLevel > 2) throw std::out_of_range("wrong CutLevel pcal FID, VW");
		
			float cutValues[3] = {17,19,22};
			float cutValuesU[3] = {390,395,400};
			
			return ((V > cutValues[cutLevel]) && (W > cutValues[cutLevel]) && (U < cutValuesU[cutLevel]));
		
		}
		

//////////////////////////////////////////////////////////////////////////////////
// Parameters for HTCCeff, will send later the code for it 
//	float HTCCEff[nHTCCXBins][nHTCCYBins];
//	readHTCCEff(HTCCEff);
    
/////////////////////////////////////////////////////////////////////////////////    
////////////// triangle Cut Params //////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    
    	float triangleCutParams[6][10][2];

        // reading parameters
        // take it from my home directory/INC/analysis2/params/*
        //this funxction reads the parameters file for the triangle cut

	
			void readTriangleCut(float cutParams[6][10][2], int isData){
			std::ifstream fData("../params/dataTriangleCut.dat");
			std::ifstream fSim("../params/simTriangleCut.dat");


//
			
			for(int s = 0; s < 6; s++){
				for(int t = 0; t < 10; t++){
					for(int l = 0; l < 2; l++){
						if(isData == 1) fData >> cutParams[s][t][l];
						if(isData == 0) fSim >> cutParams[s][t][l];
						
					}
				}
			}
		}
		

    // Function to initialize the triangle cut parameters
void initializeTriangleCut(int isData) {
    readTriangleCut(triangleCutParams, isData);
}
	  // applying the cut
	  // the input is kin.ecinE/kin.momentum, kin.pcalE/kin.momentum
	  //so ecinE is out_tree.Branch("ecinE", &ecinE) dived by electron momentum
	  // change the variable name to ecin_sf or something similar
	  // it is true for pcalE as well: out_tree.Branch("pcalE", &pcalE) divided by momentum
	  // shift is 0
	  // int pBin = (int)(p4_ele[i].P()/1); so it is binned in 1 GeV
	  
          bool SFTriangleCut(double ecinE, double pcalE, float cutParams[6][10][2], int sector, int pBin, float shift){
            sector = sector-1;
            if (pBin < 5) return true;
            //for system:
            //if (abs(shift)>0.0001) return true;
                
			float zero = cutParams[sector][pBin][0];
			float slope = cutParams[sector][pBin][1];
			if (zero + slope*ecinE + shift< pcalE) return true;
			else return false;
	  }
	  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
// helper functionc
