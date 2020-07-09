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
    string in_filename_A2 = path + "/observables/A2.txt";
    string in_filename_B2 = path + "/observables/B2.txt";
    string in_filename_A4 = path + "/observables/A4.txt";
    string in_filename_B4 = path + "/observables/B4.txt";
    ifstream in_obs_A2(in_filename_A2);
    ifstream in_obs_B2(in_filename_B2);
    ifstream in_obs_A4(in_filename_A4);
    ifstream in_obs_B4(in_filename_B4);
    
    if(!in_obs_A2 || !in_obs_B2)
    {
        cerr << "Error: file " + in_filename_A2 + " or " + in_filename_B2 + " could not be opened." << endl;
        return 1;
    }
    if(!in_obs_A4 || !in_obs_B4)
    {
        cerr << "Error: file " + in_filename_A4 + " or " + in_filename_B4 + " could not be opened." << endl;
        return 1;
    }
    
    // Open output files 
    string out_filename_A2_dim = path + "/observables/A2_dim" + cc_to_name(n) + ".txt";
    string out_filename_B2_dim = path + "/observables/B2_dim" + cc_to_name(n) + ".txt";
    string out_filename_A4_dim = path + "/observables/A4_dim" + cc_to_name(n) + ".txt";
    string out_filename_B4_dim = path + "/observables/B4_dim" + cc_to_name(n) + ".txt";
    ofstream out_obs_A2_dim(out_filename_A2_dim);
    ofstream out_obs_B2_dim(out_filename_B2_dim);
    ofstream out_obs_A4_dim(out_filename_A4_dim);
    ofstream out_obs_B4_dim(out_filename_B4_dim);

    if(!out_obs_A2_dim || !out_obs_B2_dim)
    {
        cerr << "Error: file " + out_filename_A2_dim + " or " + out_filename_B2_dim + " could not be opened." << endl;
        return 1;
    }
    if(!out_obs_A4_dim || !out_obs_B4_dim)
    {
        cerr << "Error: file " + out_filename_A4_dim + " or " + out_filename_B4_dim + " could not be opened." << endl;
        return 1;
    }

    double temp1, temp2, temp3;
    while(in_obs_A2 >> temp1 >> temp2 >> temp3)
        out_obs_A2_dim << temp1 << " " << temp2/pow(sm.dim, n) << " " << temp3/pow(sm.dim, n) << endl;
    while(in_obs_B2 >> temp1 >> temp2 >> temp3)
        out_obs_B2_dim << temp1 << " " << temp2/pow(sm.dim, n) << " " << temp3/pow(sm.dim, n) << endl;
    while(in_obs_A4 >> temp1 >> temp2 >> temp3)
        out_obs_A4_dim << temp1 << " " << temp2/pow(sm.dim, n) << " " << temp3/pow(sm.dim, n) << endl;
    while(in_obs_B4 >> temp1 >> temp2 >> temp3)
        out_obs_B4_dim << temp1 << " " << temp2/pow(sm.dim, n) << " " << temp3/pow(sm.dim, n) << endl;

    in_obs_A2.close();
    in_obs_B2.close();
    in_obs_A4.close();
    in_obs_B4.close();
    out_obs_A2_dim.close();
    out_obs_B2_dim.close();
    out_obs_A4_dim.close();
    out_obs_B4_dim.close();

    //********* END ANALYSIS **********//

    return 0;
}
