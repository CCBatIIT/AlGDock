//
//  main.cpp
//  ReferenceOBC
//
//  Created by David Minh on 10/11/16.
//  Copyright © 2016 David Minh. All rights reserved.
//

#include <iostream>
#include "ObcParameters.h"
#include "ReferenceObc.h"
//#include "ObcWrapper.h"

typedef double vector3[3];

int main(int argc, const char * argv[]) {

  int numParticles = 24;
  double charges_arr[] = {0.131300, 0.147300, 0.139400, 0.157400, 0.117000, 0.067800, 0.091200, 0.424900, 0.425600, 0.483500, 0.423600, -0.109800, -0.094800, -0.207900, -0.146800, -0.151000, 0.126300, 0.936500, -0.045900, -0.074300, -0.833000, -0.711000, -0.801600, -0.495800};
  double atomicRadii_arr[] = {0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.120000, 0.170000, 0.170000, 0.170000, 0.170000, 0.170000, 0.170000, 0.170000, 0.170000, 0.170000, 0.155000, 0.150000, 0.150000, 0.150000};
  double scaleFactors_arr[] = {0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.850000, 0.720000, 0.720000, 0.720000, 0.720000, 0.720000, 0.720000, 0.720000, 0.720000, 0.720000, 0.790000, 0.850000, 0.850000, 0.850000};
  
  std::vector<double> charges(charges_arr, charges_arr + sizeof(charges_arr) / sizeof(double));
  std::vector<double> atomicRadii(atomicRadii_arr, atomicRadii_arr + sizeof(atomicRadii_arr) / sizeof(double));
  std::vector<double> scaleFactors(scaleFactors_arr, scaleFactors_arr + sizeof(scaleFactors_arr) / sizeof(double));
  
  vector3 atomCoordinates[numParticles];
  
  atomCoordinates[0][0] = 2.805550;
  atomCoordinates[0][1] = 0.089450;
  atomCoordinates[0][2] = 1.683500;
  atomCoordinates[1][0] = 2.713690;
  atomCoordinates[1][1] = 0.503260;
  atomCoordinates[1][2] = 1.735420;
  atomCoordinates[2][0] = 2.580450;
  atomCoordinates[2][1] = 0.033700;
  atomCoordinates[2][2] = 1.759150;
  atomCoordinates[3][0] = 2.487860;
  atomCoordinates[3][1] = 0.447340;
  atomCoordinates[3][2] = 1.816930;
  atomCoordinates[4][0] = 2.946940;
  atomCoordinates[4][1] = 0.422240;
  atomCoordinates[4][2] = 1.714860;
  atomCoordinates[5][0] = 2.984070;
  atomCoordinates[5][1] = 0.255700;
  atomCoordinates[5][2] = 1.683920;
  atomCoordinates[6][0] = 2.861380;
  atomCoordinates[6][1] = 0.300150;
  atomCoordinates[6][2] = 1.449780;
  atomCoordinates[7][0] = 3.095920;
  atomCoordinates[7][1] = 0.243650;
  atomCoordinates[7][2] = 1.485840;
  atomCoordinates[8][0] = 3.085850;
  atomCoordinates[8][1] = 0.363430;
  atomCoordinates[8][2] = 1.370570;
  atomCoordinates[9][0] = 3.130310;
  atomCoordinates[9][1] = 0.402180;
  atomCoordinates[9][2] = 1.526820;
  atomCoordinates[10][0] = 2.352870;
  atomCoordinates[10][1] = 0.278670;
  atomCoordinates[10][2] = 1.858320;
  atomCoordinates[11][0] = 2.737170;
  atomCoordinates[11][1] = 0.168810;
  atomCoordinates[11][2] = 1.712970;
  atomCoordinates[12][0] = 2.683540;
  atomCoordinates[12][1] = 0.398830;
  atomCoordinates[12][2] = 1.742020;
  atomCoordinates[13][0] = 2.610740;
  atomCoordinates[13][1] = 0.137850;
  atomCoordinates[13][2] = 1.754950;
  atomCoordinates[14][0] = 2.556960;
  atomCoordinates[14][1] = 0.368490;
  atomCoordinates[14][2] = 1.788380;
  atomCoordinates[15][0] = 2.774370;
  atomCoordinates[15][1] = 0.299610;
  atomCoordinates[15][2] = 1.706940;
  atomCoordinates[16][0] = 2.522120;
  atomCoordinates[16][1] = 0.234120;
  atomCoordinates[16][2] = 1.795440;
  atomCoordinates[17][0] = 2.884280;
  atomCoordinates[17][1] = 0.507750;
  atomCoordinates[17][2] = 1.481190;
  atomCoordinates[18][0] = 2.913600;
  atomCoordinates[18][1] = 0.335030;
  atomCoordinates[18][2] = 1.658250;
  atomCoordinates[19][0] = 2.926150;
  atomCoordinates[19][1] = 0.365170;
  atomCoordinates[19][2] = 1.508870;
  atomCoordinates[20][0] = 3.065150;
  atomCoordinates[20][1] = 0.341270;
  atomCoordinates[20][2] = 1.470230;
  atomCoordinates[21][0] = 2.912950;
  atomCoordinates[21][1] = 0.548670;
  atomCoordinates[21][2] = 1.368890;
  atomCoordinates[22][0] = 2.827210;
  atomCoordinates[22][1] = 0.581440;
  atomCoordinates[22][2] = 1.566630;
  atomCoordinates[23][0] = 2.405610;
  atomCoordinates[23][1] = 0.197730;
  atomCoordinates[23][2] = 1.831700;

  ObcParameters* obcParameters = new ObcParameters(numParticles,
    ObcParameters::ObcTypeII);
  obcParameters->setAtomicRadii(atomicRadii);
  obcParameters->setScaledRadiusFactors(scaleFactors);

  obcParameters->setSolventDielectric(static_cast<double>(78.5));
  obcParameters->setSoluteDielectric(static_cast<double>(1.0));
  obcParameters->setPi4Asolv(4*M_PI*2.25936);
  obcParameters->setUseCutoff(static_cast<double>(1.5));
  
  ReferenceObc *obc = new ReferenceObc(obcParameters);
  obc->setIncludeAceApproximation(true);

  vector3 inputForces[numParticles];
  double obcEnergy;

  obcEnergy = obc->computeBornEnergyForces(atomCoordinates, charges, inputForces);

  // Energy is inconsistent with MMTK version
  std::cout << "Obc energy: " << obcEnergy << std::endl;
  for (int i = 0; i < numParticles; ++i)
    std::cout << inputForces[i][0] << ' ' << inputForces[i][1] << ' ' << inputForces[i][2] << std::endl;

  delete obcParameters;
  
  return 0;
}
