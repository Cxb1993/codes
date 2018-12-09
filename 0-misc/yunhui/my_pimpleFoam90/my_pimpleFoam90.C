/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam

Description
    Large time-step transient solver for incompressible, flow using the PIMPLE
    (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/
#include "OFstream.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
	#include <iostream>
    #include <fstream>
    using namespace std;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		double pi=3.1415926;
		double yt=3.5;
		double Ur=9.0;  //reduced velocity 
		double ZETA=0.0; //structure damping ratio ζ
		double m_star=2.0; //mass ratio
		double damp=(4*pi*ZETA)/Ur;
		double spring=(2*pi/Ur)*(2*pi/Ur);
		double inst=1/(2*m_star);
		
		double Vsolid=0.0;
		double UUY1=0.0;
		double UUY2=0.0;
		double FTY=0,YI1=0,YI2=0,YJ1=0,YJ2=0,YL1=0,YL2=0,YK1=0,YK2=0;
	
	
    Info<< "\nStarting time loop\n" << endl;
	
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
		forAll(eta,celli){
		eta[celli]=0.0 ;
		}
	forAll(eta,celli){
		double SideLength=0.05;
		double HalfDiagonal=0.866*SideLength; //HalfDiagonal==mesh.x()[celli];
		double Distance=0;
		double DistanceGrid=0;
		Distance=std::sqrt((x[celli]-7.5)*(x[celli]-7.5) +(y[celli]-yt)*(y[celli]-yt) +(z[celli]-3.5)*(z[celli]-3.5));
			if (Distance+ HalfDiagonal <= 0.5)
			{
				eta[celli]=1.0;
			}
			else if (Distance-HalfDiagonal>= 0.5)
			{
				eta[celli]=0.0;
			}
			else{
				const int rows=10;
				const int cols=10;
				const int lays=10;
				int i,j,k;
				for (i=0;i<rows;i++){
					for (j=0;j<cols;j++){
						for (k=0;k<lays;k++){
						DistanceGrid=(x[celli]-SideLength/2.0+(i+0.5)*SideLength/rows-7.5)*(x[celli]-SideLength/2+(i+0.5)*SideLength/rows-7.5)+(y[celli]-SideLength/2+(j+0.5)*SideLength/cols-yt)*(y[celli]-SideLength/2+(j+0.5)*SideLength/cols-yt)+(z[celli]-SideLength/2+(k+0.5)*SideLength/lays-3.5)*(z[celli]-SideLength/2+(k+0.5)*SideLength/lays-3.5);
							if ( DistanceGrid<=0.25){
							eta[celli]+=1/(rows*cols*lays);}
						}
					}
				}				
				} 
		}
		
		forAll(bodyforce,celli){	
				bodyforce[celli].x() = eta[celli]*((0.0-U[celli].x())/(runTime.time().deltaT().value()));
				bodyforce[celli].y() = eta[celli]*((Vsolid-U[celli].y())/(runTime.time().deltaT().value()));
				bodyforce[celli].z() = eta[celli]*((0.0-U[celli].z())/(runTime.time().deltaT().value()));
			}		
		
		
		scalar dragCoeff=0.0 ;
		forAll(bodyforce,celli){
			dragCoeff += -bodyforce[celli].x()*mesh.V()[celli]/(0.5); //1.0=velocity near boundary,ρ=1.0for noDimension
			}
		int t1 = runTime.value()-0.0+runTime.time().deltaT().value();
		if (t1 % 1 == 0 )
    		{
				ofstream myfile("Cd.dat",ios::app);
        		if(myfile.good()){
                	myfile << runTime.value() <<"\t" << dragCoeff << "\n";
        			}		
			}
		scalar LiftCoeff=0.0 ;
		forAll(bodyforce,celli){
			LiftCoeff += (-bodyforce[celli].y()*mesh.V()[celli])/(0.5); //1.0=velocity near boundary,ρ=1.0for noDimension
			}
		if (t1 % 1 == 0 )
    		{
				ofstream myfile("Cl.dat",ios::app);
        		if(myfile.good()){
                	myfile << runTime.value() <<"\t" << LiftCoeff << "\n";
        			}		
			}

//4-th Runge-Kutta Algorithm
		//1-st step
		FTY=LiftCoeff;
		YI1=Vsolid;
		YI2=inst*FTY-damp*Vsolid-spring*(yt-3.5);
		UUY1=yt+0.5*(runTime.time().deltaT().value())*YI1;
		UUY2=Vsolid+0.5*(runTime.time().deltaT().value())*YI2;
		//2-nd step
		FTY=LiftCoeff+0.5*(runTime.time().deltaT().value())*LiftCoeff;
		YJ1=Vsolid;
		YJ2=inst*FTY-damp*UUY2-spring*(UUY1-3.5);
		UUY1=yt+0.5*(runTime.time().deltaT().value())*YJ1;
		UUY2=Vsolid+0.5*(runTime.time().deltaT().value())*YJ2;
		//3-rd step
		FTY=LiftCoeff+0.5*(runTime.time().deltaT().value())*LiftCoeff;
		YK1=Vsolid;
		YK2=inst*FTY-damp*UUY2-spring*(UUY1-3.5);
		UUY1=yt+(runTime.time().deltaT().value())*YK1;
		UUY2=Vsolid+(runTime.time().deltaT().value())*YK2;
		//4-th step
		FTY=LiftCoeff+(runTime.time().deltaT().value())*LiftCoeff;
		YL1=Vsolid;
		YL2=inst*FTY-damp*UUY2-spring*(UUY1-3.5);
		//UPDATE Y AND Vs
		Vsolid=Vsolid+(YI2+2.0*YJ2+2.0*YK2+YL2)*(runTime.time().deltaT().value())/6.0;
		yt=yt+Vsolid*(runTime.time().deltaT().value());
		
		
		
			ofstream myfile("yt.dat",ios::app);
        	if(myfile.good()){
                myfile << runTime.value() <<"\t" << yt << "\n";
        		}		
			
		
		forAll(U,celli){
				U[celli].x()=U[celli].x()+bodyforce[celli].x()*(runTime.time().deltaT().value());
				U[celli].y()=U[celli].y()+bodyforce[celli].y()*(runTime.time().deltaT().value());
				U[celli].z()=U[celli].z()+bodyforce[celli].z()*(runTime.time().deltaT().value());
			}
				

		
		
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
