/*
MD.c - a simple molecular dynamics program for simulating real gas properties of Lennard-Jones particles.

 Copyright (C) 2016  Jonathan J. Foley IV, Chelsea Sweet, Oyewumi Akinfenwa

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Electronic Contact:  foleyj10@wpunj.edu
 Mail Contact:   Prof. Jonathan Foley
 Department of Chemistry, William Paterson University
 300 Pompton Road
 Wayne NJ 07470

 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

// Number of particles
int N;

//  Lennard-Jones parameters in natural units!
double sigma = 1.;
double epsilon = 1.;
double m = 1.;
double kB = 1.;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23; // m^2*kg/(s^2*K)

//  Size of box, which will be specified in natural units
double L;

//  Initial Temperature in Natural Units
double Tinit; // 2;
//  Vectors!
//
const int MAXPART = 5001;
//  Position
double r[MAXPART][3];

double rT[3][MAXPART];
//  Velocity
double v[MAXPART][3];
//  Acceleration
double a[MAXPART][3];
//  Force
double F[MAXPART][3];

// atom type
char atype[10];
//  Function prototypes
//  initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();
//  update positions and velocities using Velocity Verlet algorithm
//  print particle coordinates to file for rendering via VMD or other animation software
//  return 'instantaneous pressure'
double VelocityVerlet(double dt, int iter, FILE *fp);
//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet
void computeAccelerations();
//  Numerical Recipes function for generation gaussian distribution
double gaussdist();
//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();
//  Compute total potential energy from particle coordinates
// double Potential();
//  Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();
//  Compute total kinetic energy from particle mass and velocities
double Kinetic();

double PE;

int main()
{

    //  variable delcarations
    int i;
    double dt, Vol, Temp, Press, Pavg, Tavg, rho;
    double VolFac, TempFac, PressFac, timefac;
    double KE, /*PE,*/ mvs, gc, Z;
    char trash[10000], prefix[1000], tfn[1000], ofn[1000], afn[1000];
    FILE *infp, *tfp, *ofp, *afp;

    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  WELCOME TO WILLY P CHEM MD!\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  ENTER A TITLE FOR YOUR CALCULATION!\n");
    scanf("%s", prefix);
    strcpy(tfn, prefix);
    strcat(tfn, "_traj.xyz");
    strcpy(ofn, prefix);
    strcat(ofn, "_output.txt");
    strcpy(afn, prefix);
    strcat(afn, "_average.txt");

    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  TITLE ENTERED AS '%s'\n", prefix);
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");


    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)\n");
    printf("\n  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    scanf("%s", atype);

    if (strcmp(atype, "He") == 0)
    {

        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
    }
    else if (strcmp(atype, "Ne") == 0)
    {

        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;
    }
    else if (strcmp(atype, "Ar") == 0)
    {

        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        // strcpy(atype,"Ar");
    }
    else if (strcmp(atype, "Kr") == 0)
    {

        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
    }
    else if (strcmp(atype, "Xe") == 0)
    {

        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;
    }
    else
    {

        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        strcpy(atype, "Ar");
    }
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n                     YOU ARE SIMULATING %s GAS! \n", atype);
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n\n  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN\n");
    scanf("%lf", &Tinit);
    // Make sure temperature is a positive number!
    if (Tinit < 0.)
    {
        printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
        exit(0);
    }
    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;

    printf("\n\n  ENTER THE NUMBER DENSITY IN moles/m^3\n");
    printf("  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3\n");
    printf("  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3\n");

    scanf("%lf", &rho);

    N = 5000;
    Vol = N / (rho * NA);

    Vol /= VolFac;

    //  Limiting N to MAXPART for practical reasons
    if (N >= MAXPART)
    {

        printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", MAXPART);
        exit(0);
    }
    //  Check to see if the volume makes sense - is it too small?
    //  Remember VDW radius of the particles is 1 natural unit of length
    //  and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
    //  will be initialized with an interparticle separation equal to 2xVDW radius
    if (Vol < N)
    {

        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n", N, Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }
    // Vol = L*L*L;
    // Length of the box in natural units:
    L = pow(Vol, (1. / 3));

    //  Files that we can write different quantities to
    tfp = fopen(tfn, "w"); //  The MD trajectory, coordinates of every particle at each timestep
    ofp = fopen(ofn, "w"); //  Output of other quantities (T, P, gc, etc) at every timestep
    afp = fopen(afn, "w"); //  Average T, P, gc, etc from the simulation

    int NumTime;
    if (strcmp(atype, "He") == 0)
    {

        // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
        dt = 0.2e-14 / timefac;
        //  We will run the simulation for NumTime timesteps.
        //  The total time will be NumTime*dt in natural units
        //  And NumTime*dt multiplied by the appropriate conversion factor for time in seconds
        NumTime = 50000;
    }
    else
    {
        dt = 0.5e-14 / timefac;
        NumTime = 200;
    }

    //  Put all the atoms in simple crystal lattice and give them random velocities
    //  that corresponds to the initial temperature we have specified
    initialize();

    //  Based on their positions, calculate the ininial intermolecular forces
    //  The accellerations of each particle will be defined from the forces and their
    //  mass, and this will allow us to update their positions via Newton's law
    computeAccelerations();

    // Print number of particles to the trajectory file
    fprintf(tfp, "%i\n", N);

    //  We want to calculate the average Temperature and Pressure for the simulation
    //  The variables need to be set to zero initially
    Pavg = 0;
    Tavg = 0;

    int tenp = floor(NumTime / 10);
    fprintf(ofp, "  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");
    printf("  PERCENTAGE OF CALCULATION COMPLETE:\n  [");
    // 
    for (i = 0; i < NumTime + 1; i++)
    {

        //  This just prints updates on progress of the calculation for the users convenience
        if (i == tenp)
            printf(" 10 |");
        else if (i == 2 * tenp)
            printf(" 20 |");
        else if (i == 3 * tenp)
            printf(" 30 |");
        else if (i == 4 * tenp)
            printf(" 40 |");
        else if (i == 5 * tenp)
            printf(" 50 |");
        else if (i == 6 * tenp)
            printf(" 60 |");
        else if (i == 7 * tenp)
            printf(" 70 |");
        else if (i == 8 * tenp)
            printf(" 80 |");
        else if (i == 9 * tenp)
            printf(" 90 |");
        else if (i == 10 * tenp)
            printf(" 100 ]\n");
        fflush(stdout);

        // This updates the positions and velocities using Newton's Laws
        // Also computes the Pressure as the sum of momentum changes from wall collisions / timestep
        // which is a Kinetic Theory of gasses concept of Pressure
        Press = VelocityVerlet(dt, i + 1, tfp);
        Press *= PressFac;

        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  Now we would like to calculate somethings about the system:
        //  Instantaneous mean velocity squared, Temperature, Pressure
        //  Potential, and Kinetic Energy
        //  We would also like to use the IGL to try to see if we can extract the gas constant
        mvs = MeanSquaredVelocity();
        KE = Kinetic();
        //PE = Potential();
        //printf("---> %f\n", PE);
        
        // Temperature from Kinetic Theory
        Temp = m * mvs / (3 * kB) * TempFac;

        // Instantaneous gas constant and compressibility - not well defined because
        // pressure may be zero in some instances because there will be zero wall collisions,
        // pressure may be very high in some instances because there will be a number of collisions
        gc = NA * Press * (Vol * VolFac) / (N * Temp);
        Z = Press * (Vol * VolFac) / (N * kBSI * Temp);

        Tavg += Temp;
        Pavg += Press;

        fprintf(ofp,"  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);         
    }
    
    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    Z = Pavg*(Vol*VolFac)/(N*kBSI*Tavg);
    gc = NA*Pavg*(Vol*VolFac)/(N*Tavg);
    fprintf(afp,"  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n");
    fprintf(afp," --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n");
    fprintf(afp,"  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e         %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
    
    printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD AFTER THE SIMULATION COMPLETES\n",tfn);
    printf("\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",ofn);
    printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  \n  '%s':\n",afn);
    printf("\n  AVERAGE TEMPERATURE (K):                 %15.5f\n",Tavg);
    printf("\n  AVERAGE PRESSURE  (Pa):                  %15.5f\n",Pavg);
    printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.5f\n",gc);
    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n",100*fabs(gc-8.3144598)/8.3144598);
    printf("\n  THE COMPRESSIBILITY (unitless):          %15.5f \n",Z);
    printf("\n  TOTAL VOLUME (m^3):                      %10.5e \n",Vol*VolFac);
    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);

    fclose(tfp);
    fclose(ofp);
    fclose(afp);

    return 0;
}

void initialize()
{
    int n, p, i, j, k;
    double pos;

    n = int(ceil(pow(N, 1.0 / 3)));

    pos = L / n;

    p = 0;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                if (p < N)
                {

                    r[p][0] = (i + 0.5) * pos;
                    r[p][1] = (j + 0.5) * pos;
                    r[p][2] = (k + 0.5) * pos;

                }
                p++;
            }
        }
    }

    initializeVelocities();

}

double MeanSquaredVelocity()
{

    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;
    double v2;

    for (int i = 0; i < N; i++)
    {

        vx2 = vx2 + v[i][0] * v[i][0];
        vy2 = vy2 + v[i][1] * v[i][1];
        vz2 = vz2 + v[i][2] * v[i][2];
    }
    v2 = (vx2 + vy2 + vz2) / N;

    return v2;
}

double Kinetic()
{

    double v2, kin;

    kin = 0.;
    for (int i = 0; i < N; i++)
    {

        v2 = 0.;
        for (int j = 0; j < 3; j++)
        {
            v2 += v[i][j] * v[i][j];
        }
        kin += m * v2 / 2.;
    }

    return kin;
}


void computeAccelerations()
{
    int i, j, k;
    double f, rSqd, temp0, temp1, temp2, ri0, ri1, ri2, aux0,aux1,aux2,rSqdInv,rSqd2,rSqd4,rSqd7;

    for (i = 0; i < N; i++)
    { 
        a[i][0] = 0;
        a[i][1] = 0;
        a[i][2] = 0;
    }
    for (i = 0; i < N - 1; i++)
    {
        ri0 = r[i][0];
        ri1 = r[i][1];
        ri2 = r[i][2];
        for (j = i + 1; j < N; j++)
        {
            rSqd = 0;

            temp0 = ri0 - r[j][0];

            rSqd += temp0 * temp0;

            temp1 = ri1 - r[j][1];

            rSqd += temp1 * temp1;

            temp2 = ri2 - r[j][2];

            rSqd += temp2 * temp2;

            rSqdInv = 1.0/rSqd;
            rSqd2 = rSqdInv*rSqdInv;
            rSqd4 = rSqd2*rSqd2;
            rSqd7 =  rSqd4*rSqd2*rSqdInv;

            f = 24 * (2 * rSqd7 - rSqd4);

            aux0 = temp0 * f;
            aux1 = temp1 * f;
            aux2 = temp2 * f;
            
            a[i][0] += aux0;
            a[i][1] += aux1;
            a[i][2] += aux2;

            a[j][0] -= aux0;
            a[j][1] -= aux1;
            a[j][2] -= aux2;
            
        }
    } 
}


void computeAccelerationsPotential() {
    double var = 8 * epsilon;
    double Pot = 0.0;
    for (int i = 0; i < N; i++)
    { 
        a[i][0] = 0.0;
        a[i][1] = 0.0;
        a[i][2] = 0.0;
    }
    #pragma omp parallel
    {
        double private_a[N][3];
        for (int i = 0; i < N; i++) {
            private_a[i][0] = 0.0;
            private_a[i][1] = 0.0;
            private_a[i][2] = 0.0;
        }
        #pragma omp for schedule(dynamic) reduction(+:Pot)
        for (int i = 0; i < N - 1; i++) {
            double ri0 = r[i][0];
            double ri1 = r[i][1];
            double ri2 = r[i][2];
            for (int j = i + 1; j < N; j++) {
                double temp0 = ri0 - r[j][0];
                double temp1 = ri1 - r[j][1];
                double temp2 = ri2 - r[j][2];

                double rSqd = temp0 * temp0 + temp1 * temp1 + temp2 * temp2;
                double rSqdInv = 1.0 / rSqd;
                double rSqd2 = rSqdInv * rSqdInv;
                double rSqd3 = rSqd2 * rSqdInv;
                double rSqd4 = rSqd2 * rSqd2;
                double rSqd6 = rSqd3 * rSqd3;
                double rSqd7 = rSqd6 * rSqdInv;
                double f = 24 * (2 * rSqd7 - rSqd4);
                double aux0 = temp0 * f;
                double aux1 = temp1 * f;
                double aux2 = temp2 * f;
                private_a[i][0] += aux0;
                private_a[i][1] += aux1;
                private_a[i][2] += aux2;
                private_a[j][0] -= aux0;
                private_a[j][1] -= aux1;
                private_a[j][2] -= aux2;
                Pot += rSqd6 - rSqd3;
            }
        }
        // Atualizar a matriz 'a' fora da região paralela
        #pragma omp critical
        {
            for (int i = 0; i < N; i++) {
                a[i][0] += private_a[i][0];
                a[i][1] += private_a[i][1];
                a[i][2] += private_a[i][2];
            }
        }
    }
    PE = Pot * var;
}

//======================================================

double VelocityVerlet(double dt, int iter, FILE *fp)
{
    int i, j, k;

    double psum = 0.;
    double aux,aux2;

    aux2 = 0.5 * dt ;
    for (i = 0; i < N; i++)
    {

        aux = aux2 * a[i][0] ;
        r[i][0] += v[i][0] * dt + aux * dt;
        v[i][0] += aux;

        aux = aux2 * a[i][1]  ;
        r[i][1] += v[i][1] * dt + aux * dt;
        v[i][1] += aux;

        aux = aux2 * a[i][2]  ;
        r[i][2] += v[i][2] * dt + aux * dt;
        v[i][2] += aux;

    }
    
    computeAccelerationsPotential();
    
    for (i = 0; i < N; i++)
    {
        v[i][0] += 0.5 * a[i][0] * dt;
        v[i][1] += 0.5 * a[i][1] * dt;
        v[i][2] += 0.5 * a[i][2] * dt;
    }

    for (i = 0; i < N; i++)
    {

        if (r[i][0] < 0.)
        {
            v[i][0] *= -1.;                     
            psum += 2 * m * fabs(v[i][0]) / dt; 
        }
        if (r[i][0] >= L)
        {
            v[i][0] *= -1.;                    
            psum += 2 * m * fabs(v[i][0]) / dt;
        }

        if (r[i][1] < 0.)
        {
            v[i][1] *= -1.;                    
            psum += 2 * m * fabs(v[i][1]) / dt;
        }
        if (r[i][1] >= L)
        {
            v[i][1] *= -1.;                    
            psum += 2 * m * fabs(v[i][1]) / dt;
        }
        
        if (r[i][2] < 0.)
        {
            v[i][2] *= -1.;                    
            psum += 2 * m * fabs(v[i][2]) / dt;
        }
        if (r[i][2] >= L)
        {
            v[i][2] *= -1.;                    
            psum += 2 * m * fabs(v[i][2]) / dt;
        }

    }

    return psum / (6 * L * L);
}

void initializeVelocities()
{

    int i;

    double vCM[3] = {0, 0, 0};

    for (i = 0; i < N; i++)
    {
        v[i][0] = gaussdist();
        v[i][1] = gaussdist();
        v[i][2] = gaussdist();
        vCM[0] += m * v[i][0];
        vCM[1] += m * v[i][1];
        vCM[2] += m * v[i][2];

    }

    for (i = 0; i < 3; i++)
        vCM[i] /= N * m;

    double vSqdSum, lambda;
    vSqdSum = 0.;
    for (i = 0; i < N; i++)
    {
        v[i][0] -= vCM[0];
        v[i][1] -= vCM[1];
        v[i][2] -= vCM[2];
        vSqdSum += v[i][0] * v[i][0];
        vSqdSum += v[i][1] * v[i][1];
        vSqdSum += v[i][2] * v[i][2];

    }

    lambda = sqrt(3 * (N - 1) * Tinit / vSqdSum);

    for (i = 0; i < N; i++)
    {

        v[i][0] *= lambda;
        v[i][1] *= lambda;
        v[i][2] *= lambda;
    }
}

double gaussdist()
{
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available)
    {
        do
        {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;

        return v2 * fac;
    }
    else
    {

        available = false;
        return gset;
    }
}
