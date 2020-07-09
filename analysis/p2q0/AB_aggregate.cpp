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
    if(argc < 4)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Name of observable" << endl;
        cerr << "2) Coupling constant" << endl;
        cerr << "3) A bunch of dataset" << endl;
        return 1;
    }

    // Some declarations for later
    vector<string> datasets;
    for(int i=3; i<argc; ++i)
        datasets.push_back(argv[i]);
    double g2 = stod(argv[2]);
    string name = argv[1];
    string path = "../../data/p2q0";
    
    // Output file
    string out_filename = path + "/AB_modes/" + name + "_" + cc_to_name(g2) + ".txt";
    ofstream out_obs;
    out_obs.open(out_filename);

    if(!out_obs)
    {
        cerr << "Error: file " + out_filename + " could not be opened." << endl;
        return 1;
    }


    for(const auto& dataset : datasets)
    {

        //********* BEGIN PARAMETER INITIALIZATION **********//
        
        // Read simulation parameters from file path/dataset/init.txt
        string init_filename = path + "/" + dataset + "/init.txt";

        struct Simul_params sm;
        ifstream in_init;
        in_init.open(init_filename);

        if(!read_init_stream(in_init, sm))
        {
            cerr << "Error: couldn't read file " + init_filename << endl;
            return 1;
        }

        if(!params_validity(sm))
        {
            cerr << "Error: file " + init_filename + " does not contain the necessary parameters." << endl;
            return 1;
        }

        in_init.close();

        //********* END PARAMETER INITIALIZATION **********//

        
        //********* BEGIN READ FROM FILE **********//
        
        // Read from file path/dataset/observables/name.txt
        string in_filename = path + "/" + dataset + "/observables/" + name + ".txt";

        ifstream in_obs;
        in_obs.open(in_filename);

        if(!in_obs)
        {
            cerr << "Error: couldn't read file " + in_filename << endl;
            return 1;
        }

        double in_g2, in_val, in_err;
        while(in_obs >> in_g2 >> in_val >> in_err)
        {
            if(in_g2 == g2)
                out_obs << sm.dim << " " << in_val << " " << in_err << endl;
        }


        in_obs.close();

        //********* END READ FROM FILE **********//
    }
    
    out_obs.close();

    return 0;
}
