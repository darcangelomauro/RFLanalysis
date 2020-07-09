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
    if(argc < 5)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Name of the folder containing the data" << endl;
        cerr << "2) First index of the jobs array" << endl;
        cerr << "3) Number of jobs in the array" << endl;
        cerr << "4) Decimation factor" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    int fst_jarr = stoi(argv[2]);
    int num_jarr = stoi(argv[3]);
    int dec = stoi(argv[4]); 



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



    //********* BEGIN DECIMATION **********//
    

    // Cycle on g2 values
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Print value of g2 being precessed
        clog << "g2: " << g2 << endl;

        // Cycle on jobs in the array
        for(int i=0; i<num_jarr; ++i)
        {
            // Open input files
            string array_path = path + "/" + to_string(i+fst_jarr);
            string filename = filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            ifstream in_s, in_hl;
            string full_name_s = array_path + "/" + filename + "_S.txt";
            string full_name_hl = array_path + "/" + filename + "_HL.txt";
            in_s.open(full_name_s);
            in_hl.open(full_name_hl);
    
            // Open output file 
            ofstream out_s, out_hl;
            out_s.open(full_name_s + ".tmp");
            out_hl.open(full_name_hl + ".tmp");

            if(!(out_s && out_hl))
            {
                cerr << "Error: decimated file could not be opened." << endl;
                return 1;
            }
            
            // Create vector of correlated samples
            int length = n_meas(sm.iter_simul, sm.gap);

            // Cycle on samples
            for(int j=0; j<length; ++j) 
            {
                Geom24 G(sm.p, sm.q, sm.dim, g2);
                double S2, S4;
                in_s >> S2 >> S4;
                G.read_mat(in_hl);

                // ***** COPY DATA HERE *****
                if( !(j%dec) )
                {
                    // print S2 and S4
                    out_s << S2 << " " << S4 << endl;
                    
                    // print mat
                    for(int j=0; j<G.get_nHL(); ++j)
                    {
                        for(int k=0; k<sm.dim; ++k)
                        {
                            for(int l=0; l<sm.dim; ++l)
                                out_hl << G.get_mat(j)(k,l).real() << " " << G.get_mat(j)(k,l).imag() << " ";
                        }
                        out_hl << endl;
                    }
                }

                // ***** THAT'S IT, YOU'RE DONE *****

            }
            in_s.close();
            in_hl.close();
            out_s.close();
            out_hl.close();

            remove(full_name_s.c_str());
            rename((full_name_s+".tmp").c_str(), full_name_s.c_str());
            remove(full_name_hl.c_str());
            rename((full_name_hl+".tmp").c_str(), full_name_hl.c_str());
        }

        g2 += sm.g2_step;
    }


    
    // Log the decimation on file
    string log_filename = path + "/decimation.log";
    ofstream out_log;
    out_log.open(log_filename);

    if(!out_log)
    {
        cerr << "Error: couldn't open file " + log_filename << endl;
        return 1;
    }

    out_log << "Decimated data." << endl;
    out_log << "Decimation factor: " << dec;

    out_log.close();
    
    //********* END DECIMATION **********//

    return 0;
}
