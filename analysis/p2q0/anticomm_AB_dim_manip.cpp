#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "geometry.hpp"
#include "utils.hpp"
#include "params.hpp"
#include "statistics.hpp"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    // Check arguments
    if(argc < 2)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Path to folder containing the data" << endl;
        cerr << "2) Power of dim to divide by" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    double n = stod(argv[2]);



    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file path/init.txt
    string init_filename = path + "/init.txt";

    struct Simul_params sm;
    ifstream in_init;
    in_init.open(init_filename);

    if(!read_init_stream(in_init, sm))
    {
        cerr << "Error: couldn't read file " + init_filename << endl;
        return 1;
    }

    cout << "File " + init_filename + " contains the following parameters:" << endl;
    cout << sm.control << endl;

    if(!params_validity(sm))
    {
        cerr << "Error: file " + init_filename + " does not contain the necessary parameters." << endl;
        return 1;
    }

    in_init.close();

    //********* END PARAMETER INITIALIZATION **********//

    
    //********* BEGIN ANALYSIS **********//
    

    // Open input files
    string in_filename_A2 = path + "/observables/anticomm_AB.txt";
    ifstream in_obs_A2(in_filename_A2);
    
    if(!in_obs_A2)
    {
        cerr << "Error: file " + in_filename_A2 +" could not be opened." << endl;
        return 1;
    }
    
    // Open output files 
    string out_filename_A2_dim = path + "/observables/anticomm_AB_dim" + cc_to_name(n) + ".txt";
    ofstream out_obs_A2_dim(out_filename_A2_dim);

    if(!out_obs_A2_dim)
    {
        cerr << "Error: file " + out_filename_A2_dim + " could not be opened." << endl;
        return 1;
    }

    double temp1, temp2, temp3;
    while(in_obs_A2 >> temp1 >> temp2 >> temp3)
        out_obs_A2_dim << temp1 << " " << temp2/pow(sm.dim, n) << " " << temp3/pow(sm.dim, n) << endl;

    in_obs_A2.close();
    out_obs_A2_dim.close();

    //********* END ANALYSIS **********//

    return 0;
}
