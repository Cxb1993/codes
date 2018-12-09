#include <iostream>
#include <fstream>
#include <iomanip>          //for std::setprecision()
#include <cmath>

void meshplotter(double xLimit, double yLimit, double meshsize);
void meshgenerator(double xLimit, double yLimit, double meshsize);

int main()
{
    meshplotter(1, 1, 100);
    meshgenerator(1, 1, 100);
    return 0;
}

//Generates "xygridcoordinates.dat" file with a grid of x,y coordinates
void meshgenerator(double xLimit, double yLimit, double meshsize)
{
    using namespace std;
    
        ofstream outfile;
        outfile.open("xygridcoordinates.dat");
        
        double xInterval = xLimit/meshsize;
        double yInterval = yLimit/meshsize;
        
        for (double y = 0; y <= 1 + yInterval; y += yInterval)
        {
            for (double x = 0; x <= 1 + xInterval; x += xInterval)
                outfile << x << "," << y << "   ";
            outfile << endl;
        }
        outfile.close();
}

//Generates "xygrid.dat" file with x,y coordinates listed in column
//format for ease of plotting using gnuplot
void meshplotter(double xLimit, double yLimit, double meshsize)
{
    using namespace std;
    
        ofstream outfile;
        outfile.open("xygrid.dat");
        
        double xInterval = xLimit/meshsize;
        double yInterval = yLimit/meshsize;
        
        //Writing the data for the horizontal lines
        for (double y = 0; y < 1 + yInterval/2; y += yInterval)
        {
            for (double x = 0; x < 1 + xInterval/2; x += xInterval)
                outfile << x << " " << y << endl;
            
            //adding a linebreak to indicate next line in gnuplot
            outfile << endl;
        }

        outfile << endl << "#End of horizontal lines data" << endl <<endl;
        
        //Writing the data for the vertical lines
        for (double y = 0; y < 1 + yInterval/2; y += yInterval)
        {
            for (double x = 0; x < 1 + xInterval/2; x += xInterval)
                outfile << y << " " << x << endl;
            
            //adding a linebreak to indicate next line in gnuplot
            outfile << endl;
        }
        
        outfile << endl << "#End of horizontal lines data" << endl <<endl;
        outfile.close();
}